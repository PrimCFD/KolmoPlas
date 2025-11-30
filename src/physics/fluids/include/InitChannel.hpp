// InitChannel.hpp
#pragma once

#include "master/RunContext.hpp"
#include "master/plugin/Action.hpp"

namespace fluids
{

class InitChannel : public core::master::plugin::IAction
{
  public:
    InitChannel(double Ly, double Ubulk, double pert, void* mpi_comm);
    const core::master::plugin::ActionInfo& info() const override { return info_; }
    void execute(const core::master::MeshTileView& tile, core::master::FieldCatalog& fields,
                 double t) override;

  private:
    core::master::plugin::ActionInfo info_;
    double Ly_;
    double Ubulk_;
    double pert_;
    void* mpi_comm_{nullptr};
};

std::shared_ptr<core::master::plugin::IAction>
make_init_channel(const core::master::plugin::KV& kv, const core::master::RunContext& rc);

} // namespace fluids
