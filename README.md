# HTW_SEMESTER_PROJECT

# 🔌 Integration of Inverter-Based Resources into IEEE 13-Bus Feeder

![MATLAB/Simulink](https://img.shields.io/badge/-MATLAB%2FSimulink-0076A8?logo=mathworks&logoColor=white)
![BESS](https://img.shields.io/badge/-Battery%20Storage-4B0082)
![IEEE 13-Bus](https://img.shields.io/badge/-IEEE%2013--Bus-32CD32)

**Advanced modeling of PV systems and battery storage for grid stability in the IEEE 13-bus feeder.**  
Developed at *HTW Berlin - University of Applied Sciences* under the supervision of **Prof. Horst Schulte**.

---

## 📌 Overview
This repository contains the full report and simulation files for a semester project focused on **integrating renewable energy resources into a power distribution system**. Key components include:
- **IEEE 13-bus feeder** with unbalanced loads, voltage regulators, and mixed line configurations.
- **PV array model** with Maximum Power Point Tracking (MPPT).
- **Battery Energy Storage System (BESS)** with robust State of Charge (SOC) management.
- **Three-level NPC inverter** and harmonic filtering for grid compliance.

---

## 🚀 Key Features

### 1. **IEEE 13-Bus Feeder Implementation**
   - 4.16 kV distribution network with:
     - Overhead/underground lines
     - Delta-Star (Δ-Y) transformer (115 kV/4.16 kV)
     - Unbalanced spot/distributed loads

### 2. **PV System with MPPT**
   - Temperature/irradiance-dependent model:
     \[
     I_{pv} = [I_{pv,n} + K_I(T-T_n)]\frac{G}{G_n}
     \]
   - Open-loop MPPT for DC voltage optimization (1.2 MW capacity, 345 parallel strings).

### 3. **Battery Storage System**
   - **SOC Calculation**:
     \[
     SOC(\%) = \frac{E_{\text{kWh}}}{E_{\text{rated}}} \times 100
     \]
   - Overcharge/undercharge protection logic (90% upper limit, 10% lower limit).
   - 500 kW peak discharge | 96% round-trip efficiency.

### 4. **Grid Integration**
   - **VSC Control**: PI-based DC voltage regulation with disturbance rejection.
   - **Three-Level NPC Inverter**: <5% THD | 96% efficiency.
   - Power injection at **Bus 634 (PV)** and **Bus 635 (BESS)**.

---

## 📊 Simulation Highlights
| **Component**               | **Performance**                                                                 |
|-----------------------------|---------------------------------------------------------------------------------|
| **PV System**                | 1184 V DC output @ 1000 W/m² | 12% current rise with irradiance scaling.     |
| **BESS**                     | 80%→95% SOC in 15 mins | <10% voltage ripple | 500 kW peak discharge.          |
| **Grid Stability**           | Voltage regulation ±2% under 50% load swings | THD <3% with LC harmonic filters.            |

![System Architecture](https://via.placeholder.com/600x300?text=IEEE+13-Bus+Feeder+with+PV+and+BESS+Integration)  
*Architecture diagram from Section 5.5*

---

## 🛠️ Technical Insights

### 1. **Battery Management Logic**
```matlab
% Overcharge prevention
if SOC > 90 && P_in < 0
    P_out = 0; % Block charging
end

% Undercharge logic
if SOC < 10 && P_in > 0
    enable_charging(); % Auto-recharge at 10% SOC
end

📁 Semester_Project_HTW/
├── 📄 semester_htw_report.pdf    # Full report (Ch 1–6 + Appendix)
├── 📂 Simulink_Models/           # MATLAB models (IEEE 13-bus, PV, BESS)
├── 📂 Figures/                   # Simulation graphs (e.g., SOC, power flow)
├── 📂 Scripts/                   # Control algorithms (MPPT, VSC tuning)
└── 📂 Appendix/                  # Unit conversions & parameter tables

