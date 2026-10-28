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
- `ChActuatorSmoothing.h`: shared passive, speed-based transmission-loss crossover.
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
| `gear_loss` | Nonnegative transmission loss [W] for electric models; aggregate map conversion loss [W] for tabulated models |
| `command_saturated` | An enabled normalized command was clipped |
| `voltage_saturated` | An enabled electric voltage request exceeded the bus voltage |
| `current_saturated` | The effective current cap is active (electric model) |
| `effort_saturated` | The effort-derived current cap is active and no looser than `max_current` (electric model); both flags can be true |

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
| `max_effort` | Delivered machine-plus-transmission effort magnitude ceiling, before external friction [N or N m] | Positive infinity | Positive; infinity disables the effort-derived current cap |
| `coulomb_friction` | \f$F_c\f$ [N or N m] | 0 | Nonnegative |
| `viscous_friction` | \f$b\f$ [N s/m or N m s/rad] | 0 | Nonnegative |
| `friction_smoothing_speed` | \f$v_s\f$ [m/s or rad/s] | 0.001 | Positive |
| `allow_inconsistent_constants` | Opt out of SI constant consistency check | false | Boolean |
| `loss_smoothing_speed` | \f$w_s\f$ [m/s or rad/s], output-side transmission crossover width | 0.01 | Positive |

In coherent SI units a reciprocal motor requires \f$K_t=K_e\f$: electromagnetic power
\f$K_ti\omega_{\rm motor}\f$ must equal back-EMF power \f$K_e\omega_{\rm motor}i\f$.
This also applies to a linear motor: \f$K_t\f$ [N/A] and \f$K_e\f$ [V s/m] have the same
numerical value because force times speed equals voltage times current. By default the constructor
checks agreement to relative tolerance \f$10^{-6}\f$. `allow_inconsistent_constants=true` is an
explicit escape hatch for fitted/nonphysical data, **not** a correction: the resulting model need
not conserve energy. Convert datasheet constants first; for example, a back-EMF constant in
V/rpm is multiplied by \f$60/(2\pi)\f$ to obtain V/(rad/s), while a speed constant
\f$K_v\f$ in rpm/V gives \f$K_e=60/(2\pi K_v)\f$.

## Deriving a passive transmission loss

Let \f$e_0=NK_ti\f$ be ideal machine effort referred to the output, and
\f$P_{\rm machine}=e_0v\f$. Write delivered effort **before external friction** as
\f$e_d=e_0-\ell\f$, where \f$\ell\f$ is signed loss effort.
For motoring, efficiency means \f$e_dv=\eta_m e_0v\f$; for generating it means
\f$e_0v=\eta_g e_dv\f$ (both powers negative). Therefore the exact directional laws are
\f[
 e_d=\begin{cases}\eta_m e_0,&e_0v>0,\\ e_0/\eta_g,&e_0v<0,\end{cases}
 \qquad
 \ell=\begin{cases}(1-\eta_m)e_0,&e_0v>0,\\
                  (1-1/\eta_g)e_0,&e_0v<0.\end{cases}
\f]
Define nonnegative coefficients \f$c_m=1-\eta_m\f$ and \f$c_g=1/\eta_g-1\f$.
Then \f$\ell=c_m|e_0|\operatorname{sign}(v)\f$ in motoring and
\f$\ell=c_g|e_0|\operatorname{sign}(v)\f$ in generating.
In either case \f$P_{\rm gear}=\ell v\geq0\f$ and
\f$P_{\rm machine}=e_dv+P_{\rm gear}\f$.

At rest, efficiency ratios do not specify a unique static torque: both mechanical powers are zero.
An idealized static-friction transmission admits the interval between
\f$\eta_m e_0\f$ and \f$e_0/\eta_g\f$:
\f[
 e_d\in[\min(\eta_m e_0,e_0/\eta_g),\ \max(\eta_m e_0,e_0/\eta_g)].
\f]
**Sort the endpoints for negative effort.** A true static model would choose a value from this
band through load equilibrium or a friction constraint. The present model instead selects a
single smooth, memoryless continuation; its rest value is not a measured breakaway/stall rating.
It does not solve static-friction complementarity.

### Speed-based crossover and bounded slope

`ActuatorTransmissionLoss(e0, v, ws, cm, cg)` implements the following shared law.
Define the clamped cubic
\f[
 h(z)=\begin{cases}0,&z\leq0,\\ z^2(3-2z),&0<z<1,\\1,&z\geq1,\end{cases}
 \qquad q=|v/w_s|,
\f]
\f[
 c_d=\begin{cases}c_m,&\operatorname{sign}(e_0)=\operatorname{sign}(v),\\
                  c_g,&\text{otherwise},\end{cases}
 \qquad
 s(v)=\operatorname{sign}(v)h(2q),\qquad c(v)=c_dh(q),
\f]
\f[
 \ell=c(v)|e_0|s(v)
      =c_d|e_0|\operatorname{sign}(v)h(q)h(2q).
\f]
Set \f$s(0)=c(0)=0\f$; sign choices at zero do not affect the loss.
Both cubic factors ramp from zero immediately: there is **no central deadzone**.
The signed factor \f$s\f$ reaches \f$\operatorname{sign}(v)\f$ at \f$|v|=w_s/2\f$,
while \f$c\f$ reaches its directional coefficient at \f$|v|=w_s\f$.
At \f$|v|\geq w_s\f$, the motoring and generating laws above are **exact**, regardless of effort
magnitude. At rest \f$\ell=0\f$ and \f$e_d=e_0\f$.

Only the current quadrant's coefficient is used; there is **no averaging or cross-quadrant
coefficient contamination**, even for highly asymmetric efficiencies.
Since \f$0\leq h(q)h(2q)\leq1\f$, motoring loss satisfies
\f$|\ell|\leq c_m|e_0|\f$ and generating loss satisfies \f$|\ell|\leq c_g|e_0|\f$.
Consequently, for nonzero \f$e_0\f$, \f$e_d/e_0\in[\eta_m,1]\f$ in motoring and
\f$e_d/e_0\in[1,1/\eta_g]\f$ in generating. Thus smoothing cannot reverse motoring torque;
the signed effort lies between the corresponding ideal and far-branch values.

Since \f$c\geq0\f$ and \f$s(v)v\geq0\f$, the loss is passive everywhere, including the crossover.
For fixed \f$e_0\f$, both \f$s\f$ and \f$c\f$ are \f$C^1\f$ in speed, including zero and
the half/full-band joins. In particular, even though \f$c_d\f$ switches at a speed reversal,
\f$c=c_dh(q)\f$ rolls off quadratically to zero on both sides, with derivative zero at the crossing.
To bound the product derivative, let \f$f(q)=h(q)h(2q)\f$.
The cubic has \f$\max |h'|=3/2\f$, and the normalized signed factor has maximum slope 3.
For \f$0\leq q\leq1/2\f$, \f$h(q)\leq1/2\f$, so
\f[
 |f'(q)|\leq |h'(q)|h(2q)+h(q)\,2|h'(2q)|
           \leq\tfrac32+\tfrac12(3)=3.
\f]
For \f$q\geq1/2\f$, \f$h(2q)=1\f$ is constant and
\f$|f'(q)|=|h'(q)|\leq3/2\f$ (zero beyond \f$q=1\f$).
Thus, including the continuous derivative at zero,
\f[
 \left|\frac{\partial\ell}{\partial v}\right|
 \leq \frac{3\max(c_m,c_g)|e_0|}{w_s}.
\f]
Thus the rigorous dimensionless bound is **\f$k=3\f$**, independent of how small the effort is.
Before hard current limiting, actual delivered effort is \f$C^1\f$ in speed at fixed current,
but need not be \f$C^1\f$ in current at fixed speed (the law contains \f$|e_0|\f$).
Hard-limit transitions can also introduce derivative corners. There is no new mechanical
effort-derivative/KRM interface; the internal-current Jacobian remains the one described below.

The default `loss_smoothing_speed=1e-2` is an output-speed scale, in rad/s or m/s, **not watts**.
Choose \f$w_s\f$ **much larger than solver speed noise** and **much smaller than normal operating
speed**. When resolving a reversal, aim for \f$|\dot{v}|\Delta t\ll w_s\f$ so a step does not
skip the crossover; this is a resolution guideline, not a universal stability criterion.
Reducing it increases the worst-case velocity slope as \f$1/w_s\f$ and can worsen stiffness.
Resolve reversals and compare results after timestep refinement; do not treat smoothing as an
identified low-speed gearbox or stiction law. `friction_smoothing_speed` is a separate scale.

Output friction is
\f[
 e_f=F_c\tanh(v/v_s)+bv,\qquad e_{\rm out}=e_{\rm em}-e_f.
\f]
Friction dissipates \f$e_fv\geq0\f$ and vanishes at rest;
there is no static-friction/stiction constraint. Here \f$e_{\rm em}=e_d=e_0-\ell\f$.
A physical external Coulomb stiction model would additionally allow a rest reaction in
\f$[-F_c,F_c]\f$; the implemented \f$\tanh\f$ law instead selects zero at rest.

## Circuit, current state, and delivered-effort limits

For enabled command \f$u\f$, \f$V=V_b\operatorname{clamp}(u,-1,1)\f$; otherwise \f$V=0\f$.
Motor speed is \f$\omega_{\rm motor}=Nv\f$ and back-EMF is \f$K_eNv\f$.
With \f$L>0\f$ there is one state, \f$y_0=i\f$, initially zero, and the unconstrained circuit is
\f[
 L\dot{i}=V-Ri-K_eNv.
\f]
For current direction \f$\sigma=\pm1\f$, let
\f$\operatorname{GearLoss}(\sigma,v)\f$ denote signed loss effort evaluated with a unit current
of that sign, i.e. the helper above with \f$e_0=NK_t\sigma\f$.
Positive homogeneity in current magnitude gives the instantaneous bounds
\f[
 I_E(v,\sigma)=\frac{e_{\max}}
 {|NK_t\sigma-\operatorname{GearLoss}(\sigma,v)|},\qquad
 I_\sigma(v)=\min(I_{\max},I_E(v,\sigma)).
\f]
A zero denominator imposes no effort-derived current restriction. These bounds depend on
**speed and current direction**, not simply speed magnitude.
They enforce \f$|e_d|\leq e_{\max}\f$ without independently clipping torque.
At rest \f$I_E=e_{\max}/(NK_t)\f$. Outside the crossover, motoring has
\f$I_E=e_{\max}/(\eta_m NK_t)\f$ and generating has
\f$I_E=e_{\max}\eta_g/(NK_t)\f$.
Thus, if voltage and `max_current` permit, motoring can deliver **exactly `max_effort`**,
not the conservatively derated \f$e_{\max}\eta_m\eta_g\f$.
Braking is capped at the same delivered magnitude: `max_effort` represents a symmetric
transmission output rating, not a motor-only torque rating.

**External friction is additional** and is never clipped away: final total effort magnitude can
exceed `max_effort`. A post-hoc final torque clip could erase required braking effort while leaving
bus generation unchanged. The ODE, effort evaluation, and telemetry instead share the same current
bounds. `current_saturated` reports an active effective bound; `effort_saturated` additionally
requires the effort-derived bound to be no looser than `max_current`.

At \f$i\geq I_+(v)\f$ with positive unconstrained derivative, or
\f$i\leq-I_-(v)\f$ with negative unconstrained derivative, the RHS is zero.
An inward derivative remains active. At fixed speed, the analytic
\f$\partial\dot{i}/\partial i\f$ is \f$-R/L\f$ on the ordinary branch and zero on an
outward-clamped branch; the speed-based loss does not change it.
The clamp test uses the **strict outward unconstrained voltage residual**
\f$r=V-Ri-K_eNv\f$: \f$r>0\f$ at the upper bound or \f$r<0\f$ at the lower bound.
An equilibrium exactly on a bound has \f$r=0\f$ and retains the ordinary \f$-R/L\f$ Jacobian,
not zero merely because the RHS is zero.
Effort and telemetry **project numerical current overshoots** onto
\f$[-I_-(v),I_+(v)]\f$, but do not reset the integrated state.
Because these bounds move with speed, blocking outward derivatives is **not exact moving-bound
complementarity enforcement**: a shrinking bound may overtake the state.
Resolve crossings, commands, and speed changes with suitable timesteps and verify convergence.
Projection alone does not make an overshot state's stored energy consistent with telemetry.
This is an idealized dissipative limiter, not hardware safety protection or a detailed controller.

With \f$L=0\f$ there are no internal states:
\f[
 i_{\rm qs}=\operatorname{clamp}\left((V-K_eNv)/R,-I_-(v),I_+(v)\right).
\f]
The current dynamics are marked stiff whenever \f$L>0\f$.

**Migration:** `Parameters` and the existing model/binding names remain unchanged.
`quadrant_smoothing` has been removed entirely; replace it with `loss_smoothing_speed`.
There is no universal watts-to-speed conversion because the old transition width depended on load.
`max_effort` now rates delivered machine-plus-transmission effort **before external friction**,
using directional instantaneous current caps instead of the former conservative fixed cap.
Recheck operating points and saturation when migrating; this is not a clip of final output effort.

## Power balance and regeneration

With \f$K_t=K_e\f$ and no active limiter, the circuit and transmission balances are
\f[
 P_b=Vi=Ri^2+\frac{d}{dt}\left(\tfrac12 Li^2\right)+K_tiNv,
\f]
\f[
 P_{\rm machine}=P_m+P_{\rm gear}+P_f,\qquad
 P_{\rm gear}=\ell v,\qquad P_f=e_fv.
\f]
In motoring, energy flows bus → copper/storage → machine → transmission loss → external
friction/output. In regeneration the mechanical flow reverses, while copper, transmission, and
friction still dissipate. The full balance, including the existing ideal limiter, is
\f[
 \boxed{P_b=P_m+Ri^2+\frac{d}{dt}\left(\tfrac12 Li^2\right)
                   +P_{\rm gear}+P_f+P_{\rm lim}}.
\f]
The magnetic-storage power **\f$Li\dot{i}\f$ is essential** for transient balance; it vanishes
when \f$L=0\f$. All losses are nonnegative for a resolved admissible state with consistent constants.
The limiter term is
\f[
 P_{\rm lim}=i(V-Ri-K_eNv-L\dot{i})\geq0.
\f]
It is zero on the unconstrained circuit, and can be positive for outward-clamped dynamic current
or capped quasi-static current. **Bus telemetry remains \f$Vi\f$**: clipping current can dissipate
extra \f$P_{\rm lim}\f$; it does not silently replace bus power with mechanical power plus copper
loss. This idealized dissipation is neither `copper_loss` nor `gear_loss`, and is not a detailed
regenerative bus-controller model. Control numerical overshoots and moving-bound crossings before
using this balance quantitatively; the raw state's magnetic energy cannot be combined indiscriminately
with projected-current telemetry.
If inconsistent constants are allowed, add \f$(K_e-K_t)iNv\f$ to the right-hand side;
this term can have either sign and forfeits the physical passivity guarantee.

The four machine quadrants (nonzero current, \f$|v|\geq w_s\f$, external friction omitted) are:

| Quadrant | Speed \f$v\f$ | Current / machine effort | \f$P_{\rm machine}\f$ | Mechanical mode |
|---|---|---|---|---|
| I | + | + | + | Forward motoring |
| II | + | - | - | Forward braking/generating |
| III | - | - | + | Reverse motoring |
| IV | - | + | - | Reverse braking/generating |

Representative terminal modes, on the same far branches and without external friction, are shown
below. Torque \f$\tau=e_d\f$ has the current sign (use force for a linear actuator).

| Terminal mode | \f$V\f$ | \f$i\f$ | \f$v\f$ | \f$\tau\f$ | \f$Vi\f$ |
|---|---|---|---|---|---|
| Forward motoring | + | + | + | + | + |
| Reverse motoring | - | - | - | - | + |
| Forward regenerative braking | + | - | + | - | - |
| Reverse regenerative braking | - | + | - | + | - |
| Forward short-circuit braking | 0 | - | + | - | 0 |
| Reverse short-circuit braking | 0 | + | - | + | 0 |
| Forward plugging | - | - | + | - | + |
| Reverse plugging | + | + | - | + | + |

Regeneration requires **negative terminal power**, not merely negative output mechanical power.
In quadrant II, \f$V>0,i<0\f$ returns energy to the bus; in quadrant IV,
\f$V<0,i>0\f$ does likewise. With \f$V=0\f$, either braking quadrant is **short-circuit braking**:
bus power is zero and absorbed power becomes loss (or transient storage).
With voltage opposite speed and current in the braking direction, \f$Vi>0\f$ is **plugging**:
both the bus and mechanical system supply dissipated energy, not regeneration.
For example, positive command can regenerate when back-EMF exceeds the applied positive voltage.
Inductive transients and external friction can complicate classification; inspect signed powers
rather than inferring bus energy direction from command or final torque alone.

## Analytic checks and time scales

At rest, the steady current is
\f$i_{\rm stall}=\operatorname{clamp}(V/R,-I_-(0),I_+(0))\f$ and
\f$e_{\rm stall}=NK_ti_{\rm stall}\f$.
`GetStallEffort()` is exactly `ComputeEffort({1, true}, state_at_zero_speed)`, including the
effort-derived current bound. It is the **mathematical regularized limit as speed approaches zero
from above**, not a datasheet far-motoring-branch stall rating.
The zero-speed loss effort is zero, while a directional efficiency continuation would give
\f$\eta_m NK_ti\f$ in motoring. The discrepancy is intentional: static efficiency is indeterminate,
and this model selects \f$e_0\f$ from the static band rather than identifying a physical stall law.
Outside \f$|v|<w_s\f$, the exact motoring result is
\f$e_d=\eta_m NK_ti\f$; do not substitute it into the regularized zero-speed formula.

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
optional `double friction_smoothing_speed = 1e-3` and trailing
`double loss_smoothing_speed = 1e-2`. Existing eight- and nine-argument calls remain valid.
`SetFrictionSmoothingSpeed()` and `SetLossSmoothingSpeed()` accept only finite positive
output speeds [m/s or rad/s]; the corresponding `GetFrictionSmoothingSpeed()` and
`GetLossSmoothingSpeed()` return them. The loss scale regularizes estimated input power,
**not** the measured delivered-effort envelope.

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
    0.001, 0.01);                        // friction and conversion-loss smoothing speeds [rad/s]
~~~

If \f$E(|v|)\f$ is the interpolated envelope and \f$F(|v|)\f$ the friction curve,
\f[
 e=u_cE(|v|)-F(|v|)\tanh(v/v_s),
\f]
where \f$u_c\f$ is the clamped enabled command or zero when disabled.
The envelope limits **commanded** effort, not the total effort including friction.
Both efficiency maps are looked up at \f$(|v|,|e|)\f$, using final delivered effort.
Unlike the electric model, this effort is already an output measurement: do **not** subtract another
transmission loss from it. Instead estimate upstream input power by **adding** a passive conversion
loss, using the same speed helper but different coefficients:
\f[
 c_m^{\rm map}=1/\eta_m-1,\qquad c_g^{\rm map}=1-\eta_g,\qquad
 \ell_{\rm map}=\operatorname{ActuatorTransmissionLoss}
       (e,v,w_s,c_m^{\rm map},c_g^{\rm map}),
\f]
\f[
 P_m=ev,\qquad P_{\rm map}=\ell_{\rm map}v\geq0,\qquad P_b=P_m+P_{\rm map}.
\f]
For \f$|v|\geq w_s\f$, the map formulas are exact:
motoring gives \f$P_b=P_m+(1/\eta_m-1)P_m=P_m/\eta_m\f$;
generating gives \f$P_b=P_m-(1-\eta_g)P_m=\eta_gP_m\f$.
Inside the band, the same directional cubic-product crossover supplies a passive input-power estimate
without mixing the two map coefficients.
Mechanical delivered torque/force is unchanged at **all** speeds.
`ComputeTelemetry().efficiency` uses the signed power ratios in the telemetry section, and remains
zero at zero power or when input and output powers have opposing signs.
**Migration nuance:** tabulated `GetEfficiency()` now returns **1 at zero mechanical power**, the
limiting speed-crossover value, rather than copying the zero-power telemetry value.
Away from zero mechanical power it returns the telemetry ratio. This continuation keeps the
speed crossover continuous; it does not assert perfect physical efficiency at stall, where a power
ratio is undefined. Code needing the zero-power diagnostic convention should read telemetry.

Here `gear_loss` denotes **aggregate map conversion loss**, not an independently measured gearbox
loss. Maps must include all device losses consistently with the chosen delivered effort and friction
curve; there is no separate winding current, copper-loss estimate, or magnetic storage.
Do not add a second copper or gearbox loss to this map balance.
For fixed effort and constant map values the crossover is \f$C^1\f$ in speed, with the same slope
bound using the map coefficients. Linear/bilinear interpolation is only piecewise differentiable:
this is **not a global \f$C^1\f$ claim at table knots**, endpoint clamps, or effort-envelope corners.
In particular, disabled friction can still be reported as absorbed mechanical power with an inferred
bus return: this remains a phenomenological map, not an electrical disconnect or winding model.

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
The loss crossover uses \f$w_s=0.01\f$ rad/s. Friction is zero, the effort ceiling is infinite,
and the current bound is inactive in both runs.

The ideal no-load speed is 30 rad/s. On the exact far branches the analytic checks are
\f[
 v_m(t)=30\left(1-e^{-t/0.347222\ldots}\right),\qquad
 v_g(t)=30+10e^{-t/0.25}.
\f]
The time constants are \f$\tau_m=0.347222\ldots\f$ s in motoring and
\f$\tau_g=0.25\f$ s in generating. These solutions are not exact inside the small speed crossover,
including the initial motoring instant at rest. The initial values provide direct power checks:

| Quantity | Generating at 40 rad/s | Motoring at 0 rad/s |
|---|---|---|
| Current \f$i=(24-0.8v)/2\f$ [A] | -4 | 12 |
| Ideal output-referred machine torque \f$e_0=0.8i\f$ [N m] | -3.2 | 9.6 |
| Delivered torque \f$e_d\f$ [N m] | -4 | 9.6 |
| Bus power \f$Vi\f$ [W] | -96 | 288 |
| Copper loss \f$Ri^2\f$ [W] | 32 | 288 |
| Machine power \f$e_0v\f$ [W] | -128 | 0 |
| Transmission loss \f$\ell v\f$ [W] | 32 | 0 |
| Mechanical output \f$e_dv\f$ [W] | -160 | 0 |

Generating satisfies \f$-96=-160+32+32\f$ W: negative bus power demonstrates actual regeneration.
At rest, all mechanical powers vanish and \f$288=288\f$ W is entirely copper loss.
`GetStallEffort()` returns the regularized 9.6 N m, not the low-speed far-motoring continuation
\f$\eta_m NK_tV_b/R=8.64\f$ N m.
Outside the speed band, motoring torque is exactly \f$0.72(12-0.4v)\f$ N m;
its extrapolation toward zero speed tends to 8.64 N m, whereas the actual smooth zero-speed limit
is 9.6 N m. The two limits describe different constitutive continuations, not an inconsistency.

The console prints time [s], speed [rad/s], torque [N m], mechanical and bus power [W],
copper loss [W], current [A], and efficiency [dimensionless].
For each run it trapezoidally integrates signed bus power, copper loss, and gearbox loss, and reports
\f[
 E_{\rm residual}=\int P_b\,dt-\Delta\left(\tfrac12 Jv^2\right)
                 -\int Ri^2\,dt-\int P_{\rm gear}\,dt
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
  no-load checks; compact speed-based loss continuity, bounded derivatives, and passivity;
  asymmetric-coefficient sweeps bounding loss by its own quadrant's coefficient; current-clamp
  Jacobian and inward release; telemetry saturation flags and disabled short-circuit behavior;
  four-quadrant energy balance, speed/direction-dependent delivered-effort limiting, analytic rotor
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
