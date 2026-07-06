# Battery State-of-Charge Estimation using a Bias-Augmented Extended Kalman Filter

Estimating the State of Charge (SoC) of a lithium-ion cell in the presence of an **uncorrected current sensor bias**, comparing classical Coulomb Counting against a bias-augmented Extended Kalman Filter (EKF) implemented in MATLAB and Simulink.

## Motivation

Coulomb Counting is the simplest SoC estimation method, but it integrates current directly — any constant offset (bias) in the current sensor accumulates without bound over time. This project demonstrates that problem explicitly, and shows how augmenting an EKF's state vector with the sensor bias itself allows the filter to jointly estimate SoC **and** learn the sensor bias, using only terminal voltage measurements as feedback.

## Model

**Equivalent Circuit Model (ECM):** first-order Thevenin model — series resistance `R0` + one RC pair (`R1`, `C1`) — driven by an OCV-SoC lookup curve.

**State vector:**

```
x = [ SoC ; V1 ; bias ]
```

**Process model (discrete-time):**

```
SoC(k+1)  = SoC(k) - (eta*dt/C_nominal) * (I_meas(k) - bias(k))
V1(k+1)   = alpha*V1(k) + beta * (I_meas(k) - bias(k))
bias(k+1) = bias(k)                          (random-walk / constant model)
```

**Measurement model:**

```
Vt(k) = OCV(SoC(k)) - I(k)*R0 - V1(k)
```

The measurement Jacobian uses the local slope of the OCV-SoC curve (`dOCV/dSoC`), which is what makes the bias state observable at all — bias only affects the output indirectly, through its effect on SoC and V1.

## Repository Structure

```
battery-soc-estimation-ekf/
├── README.md
├── matlab/
│   └── battery_ekf_main.m         # Full simulation: truth model, Coulomb Counting, bias-augmented EKF
├── report/
│   └── SoC_Estimation_Report.pdf  # Short technical report with results and discussion
```

## How to Run

Requires MATLAB (no additional toolboxes beyond base MATLAB — the EKF is hand-coded).

```matlab
run('matlab/battery_ekf_main.m')
```

This will simulate 3600 s of a dynamic discharge profile, run both estimators, and produce three comparison figures (SoC estimate, SoC error, current bias estimate).


## Results

| Method             | RMSE (SoC) | Max Abs Error (SoC) |
|--------------------|:----------:|:--------------------:|
| Coulomb Counting   | 0.0344     | 0.0591               |
| Bias-Augmented EKF | 0.0066     | 0.0113               |

The EKF also recovers the injected current sensor bias (true value: 3.0 A), converging to approximately **2.5–2.7 A** within the first ~600–800 s and remaining in that band for the rest of the simulation.


## Limitations

- All results are on synthetic simulation data (no validation yet against a real battery dataset).
- The bias state is only weakly observable through the OCV-SoC slope, so the estimate does not converge exactly to the true 3.0 A value and retains some measurement-noise-driven ripple.
- The OCV-SoC curve is a simple analytic approximation, not a fitted curve from real cell characterization data.


## Author

Arooj — Electrical Engineer (Power Systems)
