Actuator Models {#manual_actuators}
================

Chrono actuator models are engine-independent components in `chrono::actuator_models`; their headers use only the C++ standard library and can be used without linking Chrono. Chrono bindings, `ChActuatorModelLinear` and `ChActuatorModelRotational`, adapt these models to `ChLinearActuator` and `ChRotationalActuator`, respectively. They are part of `Chrono_core` and require no optional module.

Positive effort extends a linear actuator or increases its rotational coordinate. Friction and other resistive terms oppose velocity. Commands are normalized; a newly constructed command is enabled by default.

Zero-state models implement `ComputeEffort` as a pure function of command and kinematic state, so repeated evaluations at the same state produce the same result. Stateful models provide initial conditions, right-hand sides, and optionally an analytic Jacobian through the same model interface. Bindings route those states through Chrono's `ChExternalDynamicsODE` integration path.

Stiff models require an implicit Chrono timestepper and a compatible direct solver. The electric actuator model includes a winding-current state when its inductance is positive; setting inductance to zero selects its quasi-static current model. Tabulated models consume user-provided speed/effort, efficiency, and friction tables and perform clamped interpolation without file I/O or bundled measurement data.
