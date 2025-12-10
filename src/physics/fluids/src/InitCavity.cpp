#include "InitCavity.hpp"

#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include "memory/MpiBox.hpp"
#include "mesh/Mesh.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

using namespace core::master;

namespace fluids
{

InitCavity::InitCavity(double Lx, double Ly, double Lz, double U0, void* mpi_comm)
    : info_{}, Lx_(Lx), Ly_(Ly), Lz_(Lz), U0_(U0), mpi_comm_(mpi_comm)
{
    info_.name = "fluids.init.cavity";
    info_.phases = plugin::Phase::PreExchange; // same as InitTG
}

std::shared_ptr<plugin::IAction> make_init_cavity(const plugin::KV& kv, const RunContext& rc)
{
    auto get = [&](const std::string& key, const std::string& def) -> std::string
    {
        auto it = kv.find(key);
        return (it == kv.end()) ? def : it->second;
    };

    auto to_d = [](const std::string& s, double fallback)
    {
        try
        {
            return std::stod(s);
        }
        catch (...)
        {
            return fallback;
        }
    };

    const double Lx = to_d(get("Lx", "1.0"), 1.0);
    const double Ly = to_d(get("Ly", "1.0"), 1.0);
    const double Lz = to_d(get("Lz", "1.0"), 1.0);
    const double U0 = to_d(get("U0", "1.0"), 1.0);

    return std::make_shared<InitCavity>(Lx, Ly, Lz, U0, rc.mpi_comm);
}

void InitCavity::execute(const MeshTileView& tile, FieldCatalog& fields, double)
{
    auto require = [&](const char* name)
    {
        if (!fields.contains(name))
        {
            throw std::runtime_error(std::string{"[fluids.init.cavity] Missing field \""} + name +
                                     "\"");
        }
    };

    require("u");
    require("v");
    require("w");
    require("p");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    if (!tile.mesh)
        throw std::runtime_error("[fluids.init.cavity] tile.mesh is null.");

    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;

    // Cell-centered (pressure) local sizes
    const int nx_c = mesh.local[0];
    const int ny_c = mesh.local[1];
    const int nz_c = mesh.local[2];

    // Staggered sizes (same convention as InitTG / VelocityCorrector)
    const int nxu = nx_c + 1;
    const int nyu = ny_c;
    const int nzu = nz_c;

    const int nxv = nx_c;
    const int nyv = ny_c + 1;
    const int nzv = nz_c;

    const int nxw = nx_c;
    const int nyw = ny_c;
    const int nzw = nz_c + 1;

    const auto ex_u = vu.extents;
    const auto ex_v = vv.extents;
    const auto ex_w = vw.extents;
    const auto ex_p = vp.extents;

    const int nxu_tot = ex_u[0];
    const int nyu_tot = ex_u[1];
    const int nzu_tot = ex_u[2];

    const int nxv_tot = ex_v[0];
    const int nyv_tot = ex_v[1];
    const int nzv_tot = ex_v[2];

    const int nxw_tot = ex_w[0];
    const int nyw_tot = ex_w[1];
    const int nzw_tot = ex_w[2];

    const int nxc_tot = ex_p[0];
    const int nyc_tot = ex_p[1];
    const int nzc_tot = ex_p[2];

    auto fail_extents = [](const char* name)
    {
        throw std::runtime_error(std::string{"[fluids.init.cavity] Unexpected extents for field "} +
                                 name);
    };

    if (nxu_tot != nxu + 2 * ng || nyu_tot != nyu + 2 * ng || nzu_tot != nzu + 2 * ng)
        fail_extents("u");
    if (nxv_tot != nxv + 2 * ng || nyv_tot != nyv + 2 * ng || nzv_tot != nzv + 2 * ng)
        fail_extents("v");
    if (nxw_tot != nxw + 2 * ng || nyw_tot != nyw + 2 * ng || nzw_tot != nzw + 2 * ng)
        fail_extents("w");
    if (nxc_tot != nx_c + 2 * ng || nyc_tot != ny_c + 2 * ng || nzc_tot != nz_c + 2 * ng)
        fail_extents("p");

    // *** THIS IS THE IMPORTANT BIT: match InitTG and use host_ptr ***
    auto* u = static_cast<double*>(vu.host_ptr);
    auto* v = static_cast<double*>(vv.host_ptr);
    auto* w = static_cast<double*>(vw.host_ptr);
    auto* p = static_cast<double*>(vp.host_ptr);

    // x-fastest indexing helpers (same layout assumptions as elsewhere)
    const auto idx_u = [=](int i, int j, int k) noexcept
    { return (k * nyu_tot + j) * nxu_tot + i; };
    const auto idx_v = [=](int i, int j, int k) noexcept
    { return (k * nyv_tot + j) * nxv_tot + i; };
    const auto idx_w = [=](int i, int j, int k) noexcept
    { return (k * nyw_tot + j) * nxw_tot + i; };
    const auto idx_p = [=](int i, int j, int k) noexcept
    { return (k * nyc_tot + j) * nxc_tot + i; };

    // Zero everything (including halos). The lid velocity will
    // be imposed by the BC handler on the first call to ApplyBCs.
    for (int k = 0; k < nzu_tot; ++k)
        for (int j = 0; j < nyu_tot; ++j)
            for (int i = 0; i < nxu_tot; ++i)
                u[idx_u(i, j, k)] = 0.0;

    for (int k = 0; k < nzv_tot; ++k)
        for (int j = 0; j < nyv_tot; ++j)
            for (int i = 0; i < nxv_tot; ++i)
                v[idx_v(i, j, k)] = 0.0;

    for (int k = 0; k < nzw_tot; ++k)
        for (int j = 0; j < nyw_tot; ++j)
            for (int i = 0; i < nxw_tot; ++i)
                w[idx_w(i, j, k)] = 0.0;

    for (int k = 0; k < nzc_tot; ++k)
        for (int j = 0; j < nyc_tot; ++j)
            for (int i = 0; i < nxc_tot; ++i)
                p[idx_p(i, j, k)] = 0.0;

    if (tile.mesh)
    {
        core::master::exchange_named_fields(fields, *tile.mesh, mpi_comm_, {"u", "v", "w", "p"});
    }
}

} // namespace fluids
