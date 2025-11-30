// InitChannel.cpp
#include "InitChannel.hpp"

#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include "memory/MpiBox.hpp"
#include "mesh/Mesh.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace core::master;

namespace fluids
{

InitChannel::InitChannel(double Ly, double Ubulk, double pert, void* mpi_comm)
    : Ly_(Ly), Ubulk_(Ubulk), pert_(pert), mpi_comm_(mpi_comm)
{
    info_.name = "init_channel";
    info_.phases =
        plugin::Phase::PreExchange; // same as InitTG :contentReference[oaicite:0]{index=0}
}

std::shared_ptr<plugin::IAction> make_init_channel(const plugin::KV& kv, const RunContext& rc)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    auto to_d = [](const std::string& s, double d)
    {
        char* e = nullptr;
        double v = std::strtod(s.c_str(), &e);
        return (e && *e == 0) ? v : d;
    };

    const double Ly = to_d(get("Ly", "1.0"), 1.0);
    const double Ubulk = to_d(get("Ubulk", "1.0"), 1.0);
    const double pert = to_d(get("pert", "0.0"), 0.0);

    return std::make_shared<InitChannel>(Ly, Ubulk, pert, rc.mpi_comm);
}

void InitChannel::execute(const MeshTileView& tile, FieldCatalog& fields, double /*t*/)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.init_channel] fields u/v/w must be registered.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");

    if (!tile.mesh)
        throw std::runtime_error("[fluids.init_channel] tile.mesh is null.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;

    const int nx_c = mesh.local[0];
    const int ny_c = mesh.local[1];
    const int nz_c = mesh.local[2];

    const int NX = mesh.global[0];
    const int NY = mesh.global[1];
    const int NZ = mesh.global[2];

    const int i0 = mesh.global_lo[0];
    const int j0 = mesh.global_lo[1];
    const int k0 = mesh.global_lo[2];

    const int nxc_tot = nx_c + 2 * ng;
    const int nyc_tot = ny_c + 2 * ng;
    const int nzc_tot = nz_c + 2 * ng;
    const int nxu_tot = (nx_c + 1) + 2 * ng;
    const int nxv_tot = nxc_tot;
    const int nyv_tot = (ny_c + 1) + 2 * ng;
    const int nxw_tot = nxc_tot;
    const int nyw_tot = nyc_tot;
    const int nzw_tot = (nz_c + 1) + 2 * ng;

    auto check = [](const char* name, const std::array<int, 3>& e, int ex, int ey, int ez)
    {
        if (e[0] != ex || e[1] != ey || e[2] != ez)
            throw std::runtime_error(std::string("[fluids.init_channel] view '") + name +
                                     "' extents do not match mesh totals.");
    };
    check("u", vu.extents, nxu_tot, nyc_tot, nzc_tot);
    check("v", vv.extents, nxv_tot, nyv_tot, nzc_tot);
    check("w", vw.extents, nxw_tot, nyw_tot, nzw_tot);

    auto idxU = [nxu_tot, nyc_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxu_tot) * (std::size_t(j) + std::size_t(nyc_tot) * std::size_t(k));
    };
    auto idxV = [nxv_tot, nyv_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxv_tot) * (std::size_t(j) + std::size_t(nyv_tot) * std::size_t(k));
    };
    auto idxW = [nxw_tot, nyw_tot](int i, int j, int k) -> std::size_t
    {
        return std::size_t(i) +
               std::size_t(nxw_tot) * (std::size_t(j) + std::size_t(nyw_tot) * std::size_t(k));
    };

    double* u = static_cast<double*>(vu.host_ptr);
    double* v = static_cast<double*>(vv.host_ptr);
    double* w = static_cast<double*>(vw.host_ptr);

    // Physical spacing in y based on global count; x/z are irrelevant for the profile
    const double dy = Ly_ / static_cast<double>(NY);

    // Laminar Poiseuille with prescribed bulk velocity:
    // U(y) = 6 * Ubulk * (y/L) * (1 - y/L)
    auto U_of_y = [&](double y)
    {
        const double eta = y / Ly_;
        return 6.0 * Ubulk_ * eta * (1.0 - eta);
    };

    // ---- U on IFaces: depends only on y (channel normal) ----
    const int nx_u = nx_c + 1;
    for (int k = 0; k < nz_c; ++k)
    {
        for (int j = 0; j < ny_c; ++j)
        {
            const double yc = (j0 + j + 0.5) * dy; // center in y
            const double U = U_of_y(yc);
            for (int i = 0; i < nx_u; ++i)
            {
                const int I = i + ng, J = j + ng, K = k + ng;
                double val = U;
                if (pert_ != 0.0)
                {
                    // cheap deterministic "noise" based on indices
                    const unsigned seed =
                        (i0 + i) * 73856093u ^ (j0 + j) * 19349663u ^ (k0 + k) * 83492791u;
                    const double r = ((seed & 0xffffu) / 32768.0) - 1.0; // [-1,1)
                    val *= (1.0 + pert_ * r);
                }
                u[idxU(I, J, K)] = val;
            }
        }
    }

    // ---- V, W = 0 everywhere ----
    for (int k = 0; k < nzc_tot; ++k)
        for (int j = 0; j < nyv_tot; ++j)
            for (int i = 0; i < nxv_tot; ++i)
                v[idxV(i, j, k)] = 0.0;

    for (int k = 0; k < nzw_tot; ++k)
        for (int j = 0; j < nyw_tot; ++j)
            for (int i = 0; i < nxw_tot; ++i)
                w[idxW(i, j, k)] = 0.0;

    // One-time halo fill
    exchange_named_fields(fields, *tile.mesh, mpi_comm_, {"u", "v", "w"});
}

} // namespace fluids
