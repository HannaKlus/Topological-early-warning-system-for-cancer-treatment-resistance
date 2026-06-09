import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d


def model(ess_vals, comp_matrix_plus, comp_matrix_P, comp_matrix_minus, time_max=20000):

    #ODE
    def model(t, y, drug_on):
        t_plus, t_P, t_minus, psa = y
        if drug_on:
            K_plus = 0.5 * t_P
            K_P = 100
        else:
            K_plus = 1.5 * t_P
            K_P = 10000
        K_minus = 10000

        r_plus = 0.0027726
        r_P = 0.0034657
        r_minus = 0.0066542

        K_plus = max(K_plus, 1e-9)


        dt_plus = r_plus * t_plus * (1 - (
                    comp_matrix_plus[0] * t_plus + comp_matrix_plus[1] * t_P + comp_matrix_plus[2] * t_minus) / K_plus)
        dt_P = r_P * t_P * (1 - (comp_matrix_P[0] * t_plus + comp_matrix_P[1] * t_P + comp_matrix_P[2] * t_minus) / K_P)
        dt_minus = r_minus * t_minus * (1 - (
                    comp_matrix_minus[0] * t_plus + comp_matrix_minus[1] * t_P + comp_matrix_minus[
                2] * t_minus) / K_minus)

        sigmaPSA = 0.5
        dpsa = (t_plus + t_P + t_minus) - sigmaPSA * psa
        return [dt_plus, dt_P, dt_minus, dpsa]


    ess_vals = np.array(ess_vals)
    ess_psa = ess_vals[3]
    PSA_zenith = ess_psa * 0.8
    PSA_nadir = PSA_zenith * 0.4


    def event_psa_drop(t, y, drug_on):
        return y[3] - PSA_nadir

    event_psa_drop.terminal = True
    event_psa_drop.direction = -1

    def event_psa_recover(t, y, drug_on):
        return y[3] - PSA_zenith

    event_psa_recover.terminal = True
    event_psa_recover.direction = 1


    initial_vals = ess_vals * 0.4
    time_current = 0
    drug_on = False

    time_points = []
    psa_data = []


    while time_current < time_max:
        current_event = event_psa_drop if drug_on else event_psa_recover

        result = solve_ivp(
            model,
            [time_current, time_max],
            initial_vals,
            args=(drug_on,),
            events=current_event,
            max_step=1,
            method="BDF"
        )

        time_points.extend(result.t)
        psa_data.extend(result.y[3])

        time_current = result.t[-1]
        initial_vals = result.y[:, -1]

        if result.status == 0:
            break

        drug_on = not drug_on


    time_points = np.array(time_points)
    psa_data = np.array(psa_data)

    time_points, unique_indices = np.unique(time_points, return_index=True)
    psa_data = psa_data[unique_indices]

    f_interp = interp1d(time_points, psa_data, kind='cubic')
    t_regular = np.arange(time_points[0], time_points[-1], 1.0)
    psa_regular = f_interp(t_regular)

    return t_regular, psa_regular