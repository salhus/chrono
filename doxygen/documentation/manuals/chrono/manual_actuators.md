Actuator Models {#manual_actuators}
================

Actuator models describe an **effort source**, not a prescribed motion or a mechanical joint.
Use them when force or torque depends on command, speed, losses, or internal actuator dynamics.
For prescribed motion and motors that include constraints, see @ref motors.

\tableofcontents

# Architecture and conventions {#actuator_models_architecture}

The headers in `src/chrono/physics/actuator_models/` implement models in
`chrono::actuator_models` using **only the C++ standard library**:

- `ChActuatorModel.h`: command, kinematic state, telemetry, and the model interface.
- `ChElectricActuatorModel.h`: a DC equivalent circuit, transmission losses, and output friction.
- `ChTabulatedActuatorModel.h`: measured effort envelopes, efficiency maps, and friction curves.

These headers require C++17 but neither Chrono headers nor a Chrono library. They perform no file I/O.
The bindings `chrono::ChActuatorModelLinear` and `chrono::ChActuatorModelRotational` belong to
`Chrono_core` (no optional module). They adapt models to `ChLinearActuator` and
`ChRotationalActuator` and integrate internal states through `ChExternalDynamicsODE`.
One immutable model can be shared by multiple bindings: each binding owns its own integrated states.

## Coordinates, signs, and units

The model sees output-side kinematics, not motor-shaft kinematics:

| Quantity | Linear binding | Rotational binding |
|---|---|---|
| `ActuatorState::position` | Extension from initial length, \f$x=s-s_0\f$ [m] | Continuous relative displacement from initial angle, \f$x=\theta-\theta_0\f$ [rad] |
| `ActuatorState::velocity` | \f$v=\dot{s}\f$ [m/s] | \f$v=\dot{\theta}\f$ [rad/s] |
| Effort \f$e\f$ | Force [N], positive in extension | Torque [N m], positive in the increasing-angle direction |
| `ActuatorState::time` | Simulation time [s] | Simulation time [s] |

`Connect()` records the reference configuration. In stand-alone/co-simulation use,
`SetActuatorInitialLength()` or `SetActuatorInitialAngle()` supplies the reference; current
coordinates are supplied with `SetActuatorLength()` or `SetActuatorAngle()`.
Supply an **unwrapped** angle in stand-alone use. Attached rotational actuators track accumulated
turns across the principal-angle boundary; updates must resolve those turns (do not skip complete
revolutions between kinematic updates). In particular, keep the magnitude of angular change between
successive updates **below \f$\pi\f$ radians** so the unwrapping direction is unambiguous.

Positive mechanical power \f$P_m=ev\f$ means delivery to the mechanical system; negative power
means mechanical absorption. Positive electrical power means drawing from the electrical bus;
negative means returning energy. Friction opposes velocity in both directions.
Neither built-in model uses `position` directly; custom models may use it for springs or travel limits.

## Commands and disabling

`ActuatorCommand` defaults to `{0, true}`. Its `value` is dimensionless, normalized to
\f$[-1,1]\f$: a fraction of bus voltage for the electric model, or a fraction of the available
effort envelope for the tabulated model. The supplied models clamp enabled commands outside this range.
Supply finite commands, kinematics, and internal states.

Bindings accept a `ChFunction` through `SetInputFunction()`; no function means **zero input**.
`GetInput(time)` clamps the function value. `SetEnabled(false)` disables the **command**, not the
physics item, and `IsEnabled()` reports that flag (initially true).

Disabling does **not** make either actuator force-free:

- Electric: terminal voltage becomes zero, but current evolves, back-EMF produces braking, and friction remains.
  This is a zero-voltage/shorted equivalent circuit, not an open-circuit disconnect.
- Tabulated: commanded effort becomes zero, but the friction curve still acts.

Do not disable the entire physics item if passive forces or internal dynamics should remain active.

# Model interface and internal states {#actuator_models_interface}

A zero-state model implements
`double ComputeEffort(const ActuatorCommand&, const ActuatorState&) const`.
It must be a **pure evaluation**: do not accumulate time, integrate state, consume a command, or
modify a hidden cache that changes results. Solvers may evaluate identical or trial states repeatedly
and in arbitrary order. A zero-state implementation returns zero from `GetNumStates()` (the default).

The optional stateful contract is:

| Method | Contract and default |
|---|---|
| `GetNumStates()` | Fixed internal-state dimension \f$n\f$; default 0. Do not change it after initialization. |
| `SetInitialConditions(y0)` | Resize/fill exactly \f$n\f$ values; default all zeros. |
| `CalculateRHS(time, y, command, state, rhs)` | Evaluate \f$\dot{y}=f(t,y,u,x,v)\f$ without mutation; resize/fill exactly \f$n\f$ derivatives. Default zeros. |
| `CalculateJacobian(time, y, command, state, J)` | Optionally fill an \f$n\times n\f$ `vector<vector<double>>`, indexed `[row][column]`, with \f$\partial f/\partial y\f$ at fixed command and mechanical state; return true. Default false requests numerical differentiation in the binding. |
| `IsStiff()` | Declare stiff internal dynamics; default false. |
| `EffortFromStates(y, command, state)` | Evaluate the output using supplied integrated states; default delegates to `ComputeEffort`. |
| `ComputeTelemetry(command, state)` | Stateless telemetry; default fills effort from `ComputeEffort` and mechanical power, leaving other fields zero/false. |
| `ComputeTelemetryFromStates(y, command, state)` | Stateful telemetry. Default preserves `ComputeTelemetry` fields and replaces effort and mechanical power using `EffortFromStates`. Override when electrical or loss quantities depend on internal states. |

The Jacobian above is **not** the complete mechanical/electrical coupled-system Jacobian.
Keep explicit `time` and `state.time` consistent in external integrations.
Bindings validate initial-condition, RHS, and successful Jacobian dimensions; invalid shapes are errors,
not requests to silently discard dynamics.

Call `Connect()` (which initializes) for attached use, or `Initialize()` for stand-alone use, before
evaluating a stateful binding. Chrono advances the binding's state vector; model evaluators do not.
Calling the electric model's `ComputeEffort` directly always evaluates its **quasi-static** current
solution, even when its configured state dimension is one. To use a supplied current, call
`EffortFromStates` and `ComputeTelemetryFromStates`.

## Telemetry

Both bindings expose `GetTelemetry(double time) const`. It evaluates the current kinematics and
command with the **current integrated states** and does not advance the system. The time argument
selects the input-function evaluation time; it does not reconstruct historical states.

| Field | Meaning |
|---|---|
| `applied_effort` | Final signed output force [N] or torque [N m], including friction and model constraints |
| `mechanical_power` | Final applied effort times output velocity [W] |
| `electrical_power` | Signed terminal/bus power [W]; electric model uses \f$Vi\f$ |
| `efficiency` | Instantaneous conversion ratio in [0,1], or the tabulated estimate |
| `current` | Electric winding current [A]; zero if not supplied by the model |
| `copper_loss` | Electric winding loss \f$Ri^2\f$ [W]; zero if not supplied |
| `command_saturated` | An enabled normalized command was clipped |
| `voltage_saturated` | An enabled electric voltage request exceeded the bus voltage |
| `current_saturated` | The effective current cap is active (electric model) |
| `effort_saturated` | The tighter effort-derived current cap is active (electric model); both flags can be true |

For the electric model the final signed conversion ratio is
\f[
 \eta_{\rm reported} =
 \begin{cases}
  \operatorname{clamp}(P_m/P_b,0,1), & P_m>0,\ P_b>0,\\
  \operatorname{clamp}(P_b/P_m,0,1), & P_m<0,\ P_b<0,\\
  0, & \text{otherwise}.
 \end{cases}
\f]
Zero powers and opposing signs report zero, not an absolute-value ratio.
Clamping a ratio does not prove passivity: stored magnetic energy can temporarily produce mechanical
output while bus power is zero or of the opposite sign. Use the energy balance below.
In a binding the input function is already clamped by `GetInput`, so model telemetry cannot recover
the original out-of-range request; inspect the original function if that distinction matters.

# Electric actuator {#actuator_models_electric}

`ChElectricActuatorModel::Parameters` uses motor-side electrical constants and output-side friction.
For rotational use let \f$N=\omega_{\rm motor}/v>0\f$. For a rotary-to-linear transmission,
\f$N=\omega_{\rm motor}/v\f$ instead has units rad/m; use the corresponding force and linear-speed units.
For a direct linear-motor equivalent, use motor-side linear speed instead, \f$K_t\f$ [N/A],
\f$K_e\f$ [V/(m/s)], and a dimensionless motor/output linear-speed ratio.
The transmission is ideal apart from the specified efficiencies; it adds no rotor inertia, backlash,
elasticity, thermal state, battery, switching converter, or bus-energy storage.

## Parameters

Unless stated otherwise, numeric parameters must be finite. Invalid parameters throw `std::invalid_argument`.

| Member | Symbol / units | Default | Allowed range |
|---|---|---|---|
| `resistance` | \f$R\f$ [ohm] | 1 | \f$R>0\f$ |
| `inductance` | \f$L\f$ [H] | 0.01 | \f$L\geq0\f$; zero selects quasi-static current |
| `torque_constant` | \f$K_t\f$ [N m/A], or [N/A] for a linear motor | 1 | Positive |
| `back_emf_constant` | \f$K_e\f$ [V/(rad/s)], or [V/(m/s)] for a linear motor | 1 | Positive |
| `bus_voltage` | \f$V_b\f$ [V] | 24 | Positive |
| `gear_ratio` | \f$N\f$ [dimensionless same-kind speed ratio, or rad/m for rotary-to-linear transmission] | 1 | Positive |
| `motoring_efficiency` | \f$\eta_m\f$ [dimensionless] | 0.9 | \f$0<\eta_m\leq1\f$ |
| `generating_efficiency` | \f$\eta_g\f$ [dimensionless] | 0.8 | \f$0<\eta_g\leq1\f$ |
| `max_current` | \f$I_{\max}\f$ [A] | 10 | Positive |
| `max_effort` | Electromagnetic output effort ceiling [N or N m] | Positive infinity | Positive; infinity disables the effort-derived current cap |
| `coulomb_friction` | \f$F_c\f$ [N or N m] | 0 | Nonnegative |
| `viscous_friction` | \f$b\f$ [N s/m or N m s/rad] | 0 | Nonnegative |
| `friction_smoothing_speed` | \f$v_s\f$ [m/s or rad/s] | 0.001 | Positive |
| `allow_inconsistent_constants` | Opt out of SI constant consistency check | false | Boolean |
| `quadrant_smoothing` | \f$\epsilon_p\f$ [W] | 0.001 | Positive |

In coherent SI units a reciprocal motor requires \f$K_t=K_e\f$. By default the constructor checks
agreement to relative tolerance \f$10^{-6}\f$. `allow_inconsistent_constants=true` is an explicit
escape hatch for fitted/nonphysical data, **not** a correction: the resulting model need not conserve
energy. Convert rpm-based constants before constructing the model.

## Circuit, current state, and smooth transmission

For enabled command \f$u\f$, \f$V=V_b\operatorname{clamp}(u,-1,1)\f$; otherwise \f$V=0\f$.
Motor speed is \f$\omega_{\rm motor}=Nv\f$ and back-EMF is \f$K_eNv\f$.
With \f$L>0\f$ there is one state, \f$y_0=i\f$, initially zero:
\f[
 L\dot{i}=V-Ri-K_eNv.
\f]
Both current and electromagnetic effort limits are enforced through the fixed bound
\f[
 I_{\rm lim}=\min\left(I_{\max},\frac{e_{\max}\eta_g}{NK_t}\right).
\f]
Since \f$1/\eta_g\f$ is the largest transmission factor, this conservative current cap guarantees
\f$|e_{\rm em}|\leq e_{\max}\f$ in every quadrant without independently clipping output torque.
For finite `max_effort`, the bound can reduce available motoring and stall effort below that ceiling.
When the effort-derived current cap dominates, achievable zero-speed stall effort is at most
\f$e_{\max}\eta_g\f$, and far-branch motoring effort is at most
\f$e_{\max}\eta_m\eta_g\f$; voltage/resistance can reduce both further.
For example, `max_effort=100`, \f$\eta_m=0.9\f$, and \f$\eta_g=0.8\f$ give ceilings of
80 at stall and 72 in far-branch motoring, in the chosen effort units, not a 100-unit stall rating.
**Passive friction is additional** and is never clipped away: the final total effort magnitude can
exceed `max_effort`. This avoids erasing dissipative back-driving torque or decoupling mechanical
effort from the electrical current.
The same effective current cap is used by the ODE, effort evaluation, and telemetry.
`effort_saturated` is true when the effort-derived cap is tighter than `max_current` and the effective
cap is active; `current_saturated` reports the effective cap regardless of which parameter sets it.

**Migration note:** `max_effort` is an electromagnetic effort ceiling **before passive friction**,
not a hard clip of final force/torque. A final torque clip can erase required braking effort while
leaving generated electrical power unchanged, creating energy during regeneration.
This conservative constitutive cap is not a hardware safety limiter. Model safety protection,
controller switching, or mechanical torque-limiting devices separately in the plant if required.

At \f$i\geq I_{\rm lim}\f$ with positive unconstrained derivative, or
\f$i\leq-I_{\rm lim}\f$ with negative unconstrained derivative, the RHS is set to zero.
An inward derivative remains active, allowing recovery from saturation.
The analytic current Jacobian is \f$-R/L\f$ normally and **exactly zero** in an outward-clamped row.
Effort and telemetry use current subject to the configured limits; a derivative clamp is not an
exact projection of arbitrary numerical overshoots, so resolve limit crossings with appropriate steps.

With \f$L=0\f$ there are no internal states and current is obtained algebraically:
\f[
 i_{\rm qs}=\operatorname{clamp}\left((V-K_eNv)/R,-I_{\rm lim},I_{\rm lim}\right).
\f]
The current dynamics are marked stiff whenever \f$L>0\f$.

To avoid a discontinuous gearbox torque at a power-quadrant transition, define
\f[
 p=NK_tiv,\qquad
 q=\min(|p|/\epsilon_p,1),\qquad h=q^2(3-2q),
\f]
\f[
 a(p)=
 \begin{cases}
  1+(\eta_m-1)h,&p\geq0,\\
  1+(1/\eta_g-1)h,&p<0.
 \end{cases}
\f]
The electromagnetic output effort before friction is \f$e_{\rm em}=NK_ti\,a(p)\f$.
This compact \f$C^1\f$ blend reaches the **exact** motoring branch \f$a=\eta_m\f$ for
\f$p\geq\epsilon_p\f$ and generating branch \f$a=1/\eta_g\f$ for
\f$p\leq-\epsilon_p\f$. At zero power \f$a=1\f$.
Thus \f$(1-a)p\geq0\f$: smoothing itself cannot create energy in either quadrant.
It is a regularization near zero power, not a measured low-speed gearbox law.

Output friction is
\f[
 e_f=F_c\tanh(v/v_s)+bv,\qquad e_{\rm out}=e_{\rm em}-e_f.
\f]
Friction dissipates \f$e_fv\geq0\f$ and vanishes at rest;
there is no static-friction/stiction constraint.

## Power balance and regeneration

Without active current/output limits, and with \f$K_t=K_e\f$,
\f[
 P_b=Vi=Ri^2+\frac{d}{dt}\left(\tfrac12 Li^2\right)+K_tiNv,
\f]
\f[
 P_b=P_m+Ri^2+Li\dot{i}+(1-a)p+e_fv.
\f]
The magnetic-storage power **\f$Li\dot{i}\f$ is essential** for transient balance.
For \f$L=0\f$ it vanishes. The extra term for inconsistent constants is
\f$(K_e-K_t)iNv\f$, which can have either sign.
The hard current constraint represents an ideal dissipative limiter, not a detailed bus controller.
For a bounded state its extra power is
\f[
 P_{\rm lim}=i(V-Ri-K_eNv-L\dot{i})\geq0.
\f]
Add this term to the right-hand side of the balance when the limiter is active (also for capped
quasi-static current). It is not included in `copper_loss`, and numerical current overshoots must be
controlled before using this balance quantitatively.

Regeneration requires **negative terminal power**, not merely negative output mechanical power.
For example, zero applied voltage gives zero bus power even while back-driving produces current,
copper loss, and mechanical braking. Conversely a positive-voltage command can regenerate when
back-EMF exceeds that voltage and current becomes negative.

## Analytic checks and time scales

At rest, the steady current is
\f$i_{\rm stall}=\operatorname{clamp}(V/R,-I_{\rm lim},I_{\rm lim})\f$ and the smooth model gives
\f$e_{\rm stall}=NK_ti_{\rm stall}\f$ (unity transmission factor at zero power).
This is the actual voltage-dependent steady stall result, not a current-limit rating.
`GetStallEffort()` evaluates this expression at the full positive bus voltage, including the
effort-derived current bound.
The absence of \f$\eta_m\f$ at zero speed is **intentional**: the smooth coefficient is unity at
zero electromagnetic power, rather than extending the far-branch efficiency to the crossover.
Away from the smoothing region, the usual motoring approximation is
\f$e\simeq\eta_mNK_ti\f$; do not substitute that approximation into the exact zero-speed formula.

The ideal full-voltage, frictionless no-load speed is
\f$v_0=V_b/(K_eN)\f$, returned by `GetNoLoadSpeed()`.
Friction, external load, and finite command change the operating speed.
At fixed speed with inactive limits,
\f[
 \tau_e=L/R,\quad
 i_\infty=(V-K_eNv)/R,\quad
 i(t)=i_\infty+[i(0)-i_\infty]\exp(-t/\tau_e).
\f]
For an unclipped, frictionless quasi-static motor driving output inertia \f$J\f$ on a **fixed far
branch** with factor \f$a=\eta_m\f$ or \f$1/\eta_g\f$,
\f[
 J\dot{v}=\frac{aNK_t}{R}(V-K_eNv),\qquad
 \tau_m=\frac{JR}{aN^2K_tK_e}.
\f]
For linear use replace \f$J\f$ by output mass. These expressions are local checks, not a global
solution across reversals, the smoothed region, saturation, friction, or varying load.
For finite inductance the coupled small-signal characteristic equation is
\f$LJ\lambda^2+RJ\lambda+aN^2K_tK_e=0\f$ under the same fixed-branch assumptions.

# Tabulated measured actuator {#actuator_models_tabulated}

`ChTabulatedActuatorModel` has zero states. Its constructor accepts eight vectors, followed by
an optional ninth argument `double friction_smoothing_speed = 1e-3`.
Existing eight-argument calls remain valid. `SetFrictionSmoothingSpeed()` accepts only finite positive
speeds [m/s or rad/s], and `GetFrictionSmoothingSpeed()` returns the value.

All data are **output-side SI quantities**; do not apply a second gear reduction:

| Argument | Units and requirements |
|---|---|
| `speed_axis`, `max_effort` | Nonnegative speed magnitudes [m/s or rad/s] and available effort magnitudes [N or N m]; equal nonzero lengths |
| `efficiency_speed_axis` | Nonnegative speed magnitudes [m/s or rad/s] |
| `efficiency_effort_axis` | Nonnegative effort magnitudes [N or N m] |
| `motoring_efficiency`, `generating_efficiency` | Dimensionless values in (0,1]; each has speed-axis length times effort-axis length entries |
| `friction_speed_axis`, `friction_effort` | Nonnegative speed and friction-effort magnitudes, with equal nonzero lengths |

Every axis must be finite, nonempty, and strictly increasing (one-point axes are valid).
All table values must be finite; effort magnitudes must be nonnegative. Invalid data throw
`std::invalid_argument`. One-dimensional curves use linear interpolation and efficiency maps use
bilinear interpolation. Queries outside an axis clamp to its endpoint; they do not extrapolate.

Efficiency entries are **row-major with speed first**:
`values[speed_index * efficiency_effort_axis.size() + effort_index]`.
For example, the following illustrative data could be replaced by measurements; Chrono supplies no
measurement files and does not parse CSV:

~~~{.cpp}
#include "chrono/physics/actuator_models/ChTabulatedActuatorModel.h"

chrono::actuator_models::ChTabulatedActuatorModel measured(
    {0.0, 10.0}, {100.0, 60.0},          // speed [rad/s], effort envelope [N m]
    {0.0, 10.0}, {0.0, 100.0},          // efficiency speed and effort axes
    {0.80, 0.85, 0.86, 0.90},           // speed=0 row, then speed=10 row
    {0.70, 0.75, 0.76, 0.80},           // generating rows, same layout
    {0.0, 10.0}, {1.0, 2.0},            // friction speed and magnitude
    0.01);                              // smoothing speed [rad/s]
~~~

If \f$E(|v|)\f$ is the interpolated envelope and \f$F(|v|)\f$ the friction curve,
\f[
 e=u_cE(|v|)-F(|v|)\tanh(v/v_s),
\f]
where \f$u_c\f$ is the clamped enabled command or zero when disabled.
The envelope limits **commanded** effort, not the total effort including friction.
Efficiency is looked up at \f$(|v|,|e|)\f$ using the sign of final \f$ev\f$ to select the map.
The inferred electrical power is \f$P_b=P_m/\eta_m\f$ when motoring and
\f$P_b=\eta_gP_m\f$ when generating; reported efficiency is zero at zero power.
There is no winding current, copper loss, or magnetic storage in this model.
The maps should describe the complete measured output/bus conversion consistently with the chosen
friction curve. In particular, disabled friction can be reported as absorbed mechanical power with an
inferred bus return: this is a phenomenological map, not an electrical disconnect or winding model.

# Using the models {#actuator_models_usage}

## External program: no Chrono library

Users who need only a constitutive actuator can include the headers directly, for example from an
installed include directory or the source tree. This complete program uses only the standard library:

~~~{.cpp}
#include <iostream>
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

int main() {
    using namespace chrono::actuator_models;
    ChElectricActuatorModel::Parameters p;
    p.inductance = 0;                    // algebraic current, no external ODE required
    ChElectricActuatorModel model(p);
    ActuatorCommand command{0.5, true};
    ActuatorState state{0.0, 2.0, 0.0};  // displacement [rad], velocity [rad/s], time [s]
    const auto t = model.ComputeTelemetry(command, state);
    std::cout << t.applied_effort << " N m, " << t.current << " A, "
              << t.electrical_power << " W bus, " << t.copper_loss << " W copper\n";
}
~~~

For example: `c++ -std=c++17 -I/path/to/chrono/src example.cpp -o example`.
For \f$L>0\f$, the external application must initialize and integrate `y` using the stateful contract;
passing time to `ComputeEffort` does not integrate current.

## Linear Chrono binding

The following complete minimal program applies a measured force to a free body along a nonzero-length
line. It intentionally does not constrain other degrees of freedom.

~~~{.cpp}
#include <memory>
#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelLinear.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChTabulatedActuatorModel.h"

int main() {
    using namespace chrono;
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    auto ground = std::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);
    auto slider = std::make_shared<ChBody>();
    slider->SetMass(1);
    slider->SetInertiaXX(ChVector3d(1, 1, 1));
    slider->SetPos(ChVector3d(1, 0, 0));
    system.AddBody(slider);

    auto model = std::make_shared<actuator_models::ChTabulatedActuatorModel>(
        std::vector<double>{0, 1}, std::vector<double>{10, 5},
        std::vector<double>{0, 1}, std::vector<double>{0, 10},
        std::vector<double>{0.9, 0.9, 0.9, 0.9},
        std::vector<double>{0.8, 0.8, 0.8, 0.8},
        std::vector<double>{0, 1}, std::vector<double>{0.1, 0.2});
    auto actuator = std::make_shared<ChActuatorModelLinear>(model);
    actuator->SetInputFunction(std::make_shared<ChFunctionConst>(0.5));
    actuator->Connect(ground, slider, true, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
    system.Add(actuator);
    system.DoStepDynamics(1e-3);
    const auto telemetry = actuator->GetTelemetry(system.GetChTime());
    return telemetry.applied_effort > 0 ? 0 : 1;
}
~~~

`local=true` makes the connection points body-local. Positive force pushes these endpoints apart.
Avoid coincident endpoints, for which the force direction is undefined.
Add a prismatic joint separately if the application requires a guided slider.

## Rotational binding with current dynamics

This complete example applies torque about the absolute actuator frame's Z axis.
The bodies are free except for the fixed ground; add a revolute joint separately for a supported shaft.

~~~{.cpp}
#include <memory>
#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

int main() {
    using namespace chrono;
    ChSystemNSC system;
    system.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    auto ground = std::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);
    auto rotor = std::make_shared<ChBody>();
    rotor->SetMass(1);
    rotor->SetInertiaXX(ChVector3d(1, 1, 1));
    system.AddBody(rotor);

    actuator_models::ChElectricActuatorModel::Parameters p;
    p.inductance = 0.01;
    auto model = std::make_shared<actuator_models::ChElectricActuatorModel>(p);
    auto actuator = std::make_shared<ChActuatorModelRotational>(model);
    actuator->SetInputFunction(std::make_shared<ChFunctionConst>(0.25));
    actuator->Connect(ground, rotor, ChFrame<>(ChVector3d(0, 0, 0)));
    system.Add(actuator);
    for (int step = 0; step < 100; ++step)
        system.DoStepDynamics(1e-4);
    const auto telemetry = actuator->GetTelemetry(system.GetChTime());
    actuator->SetEnabled(false);         // future steps retain passive behavior and current decay
    return telemetry.current > 0 ? 0 : 1;
}
~~~

For stiff internal states, select an **implicit timestepper and a compatible direct solver** before
stepping; the example uses the core sparse QR solver without an optional solver module.
An explicit scheme generally needs a step much smaller than \f$L/R\f$ (the scalar forward-Euler
stability condition is \f$\Delta t<2L/R\f$, not an accuracy recommendation).
Implicit stability does not resolve a transient automatically: start with a step well below
\f$L/R\f$ (for example \f$(L/R)/10\f$), also resolve mechanical motion, commands, friction smoothing,
and limit crossings, then halve the step to check convergence.
Choose \f$L=0\f$ only when the electrical transient is intentionally neglected; it changes the model.

Users integrating these bindings into an **external project** should follow `template_project/`,
build/install Chrono, and use `find_package(Chrono CONFIG REQUIRED)` and `${CHRONO_TARGETS}`.
`Chrono_DIR` may point to a build or install tree; a local installation is usually cleanest.
Do not copy the core binding sources into an external application.

## Headless motoring and regeneration demo

`src/demos/mbs/demo_MBS_actuator_models.cpp` runs without visualization, optional modules, input
files, or CSV-output configuration. It performs two **independent** full-voltage runs, each for 2 s
with a \f$10^{-4}\f$ s step: generating from 40 rad/s and motoring from rest.
The parameters are \f$R=2\f$ ohm, \f$L=0\f$ H, \f$K_t=K_e=0.2\f$ in SI,
\f$V_b=24\f$ V, \f$N=4\f$, \f$J=0.1\f$ kg m\f$^2\f$,
\f$\eta_m=0.9\f$, \f$\eta_g=0.8\f$, and \f$I_{\max}=20\f$ A.
Friction is zero, the effort ceiling is infinite, and the current bound is inactive in both runs.

The ideal no-load speed is 30 rad/s. On the exact far branches the analytic checks are
\f[
 v_m(t)=30\left(1-e^{-t/0.347222\ldots}\right),\qquad
 v_g(t)=30+10e^{-t/0.25}.
\f]
These are not exact inside the small smooth crossover region, including the initial motoring instant
at rest. For example, `GetStallEffort()` gives \f$4(0.2)(24/2)=9.6\f$ N m there, not the
far-motoring approximation \f$0.9(9.6)=8.64\f$ N m.
Generating starts with \f$i=(24-0.2(4)(40))/2=-4\f$ A and \f$P_b=-96\f$ W:
negative bus power demonstrates actual regeneration.

The console prints time [s], speed [rad/s], torque [N m], mechanical and bus power [W],
copper loss [W], current [A], and efficiency [dimensionless].
For each run it trapezoidally integrates signed bus power, copper loss, and gearbox loss, and reports
\f[
 E_{\rm residual}=\int P_b\,dt-\Delta\left(\tfrac12 Jv^2\right)
                 -\int Ri^2\,dt-\int(1-a)p\,dt
\f]
in joules. Since \f$L=0\f$, there is no magnetic-storage term. The residual is a numerical diagnostic,
not an exactly zero identity for a finite timestep; reduce the step to check convergence.

Developers can build and run it with:

~~~{.sh}
cmake -S . -B build -DBUILD_DEMOS=ON
cmake --build build --target demo_MBS_actuator_models --parallel 2
./build/bin/demo_MBS_actuator_models
~~~

# Developer verification and extension {#actuator_models_development}

The relevant existing test programs are:

- `utest_CH_actuator_models`: SI-constant and parameter validation; analytic current, stall, and
  no-load checks; compact quadrant smoothing continuity, derivatives, and passivity; current-clamp
  Jacobian and inward release; telemetry saturation flags and disabled short-circuit behavior;
  four-quadrant energy balance, conservative effort-derived current limiting, analytic rotor
  integration and integrated energy; table interpolation, clamping, validation, and purity;
  binding force/torque directions and implicit current integration.
- `utest_CH_actuator_model_bindings`: zero/reference positions; default-zero, clamped, null-function,
  and disabled commands through effort, RHS, Jacobian, and telemetry; empty Jacobian resizing and
  malformed-dimension rejection; connected linear displacement/current RHS; attached accumulated
  revolutions in both directions; telemetry from integrated rather than quasi-static current in
  both linear and rotational bindings.
- `utest_CH_actuator_models_isolated`: instantiates and executes **both** electric and tabulated
  models, including telemetry and disabled passive behavior, without Chrono or Eigen include
  paths/linkage.
- `utest_CH_actuators`: the underlying actuator integration path.

For **developers modifying this repository**, initialize bundled test dependencies and configure
an out-of-source build explicitly:

~~~{.sh}
git submodule init
git submodule update
cmake -S . -B build -DBUILD_TESTING=ON -DBUILD_TESTING_BASE=ON -DBUILD_DEMOS=OFF
cmake --build build --target utest_CH_actuator_models utest_CH_actuator_model_bindings utest_CH_actuator_models_isolated utest_CH_actuators --parallel 2
ctest --test-dir build -R actuator --output-on-failure
~~~

Add `--config Release` to build and CTest commands for a multi-configuration generator.
`BUILD_TESTING` defaults **OFF** at the top level; `BUILD_TESTING_BASE` defaults **ON** only within
the unit-test configuration and is an advanced cached option. Both must be enabled to generate the
physics test targets. A stale cache with either option disabled explains an "unknown target"/"no rule
to make target" error: rerun the configuration command against the **same build directory**.
In addition, `src/CMakeLists.txt` forcibly caches `BUILD_TESTING=OFF` if the bundled GoogleTest
`CMakeLists.txt` is missing. After initializing submodules, explicitly set `BUILD_TESTING=ON`
again; merely rebuilding will not restore the targets. Initializing the bundled GoogleTest/Google
Benchmark submodules is required for the corresponding test/benchmark configuration; do not replace
them with an unrelated system test framework.

When adding a model:

1. Keep its header under `src/chrono/physics/actuator_models/` standard-library-only; do not introduce
   Eigen, Chrono types, global state, file parsing, or mandatory external dependencies.
2. Specify coordinate references, units, signs, normalization, disabled/passive behavior, and all
   parameter ranges. Reject invalid parameters and tables.
3. Keep zero-state evaluations pure; place every evolving quantity in the declared state vector.
4. Implement exact state/RHS dimensions, effort from states, stiffness, and an analytic Jacobian
   where useful. Check it against finite differences, including active constraint branches.
5. Override stateful telemetry when needed; verify mechanical effort matches telemetry and account
   for losses and stored energy separately in all four effort/velocity quadrants.
6. Add focused regression cases to the existing test programs, including default-zero commands,
   disabling, saturation, zero crossings, and external header isolation. Register any new target in
   the owning `CMakeLists.txt`.
7. Verify attached and stand-alone bindings, repeated evaluations, and timestep convergence. Update
   this manual and the changelog when the public contract changes.
