# from woma
G = 6.67408e-11  # m^3 kg^-1 s^-2
R_earth = 6.371e6  # m
M_earth = 5.9724e24  # kg
R_gas = 8.3145  # Gas constant (J K^-1 mol^-1)
L_em = 3.5e34

# Material IDs, same as SWIFT ( = type_id * type_factor + unit_id)
type_factor = 100
Di_mat_type = {
    "idg": 0,
    "Til": 1,
    "HM80": 2,
    "SESAME": 3,
    "ANEOS": 4,
}
Di_mat_id = {
    # Ideal Gas
    "idg_HHe": Di_mat_type["idg"] * type_factor,
    "idg_N2": Di_mat_type["idg"] * type_factor + 1,
    "idg_CO2": Di_mat_type["idg"] * type_factor + 2,
    # Tillotson
    "Til_iron": Di_mat_type["Til"] * type_factor,
    "Til_granite": Di_mat_type["Til"] * type_factor + 1,
    "Til_water": Di_mat_type["Til"] * type_factor + 2,
    "Til_basalt": Di_mat_type["Til"] * type_factor + 3,
    # Hubbard & MacFarlane (1980) Uranus/Neptune
    "HM80_HHe": Di_mat_type["HM80"] * type_factor,  # Hydrogen-helium atmosphere
    "HM80_ice": Di_mat_type["HM80"] * type_factor + 1,  # H20-CH4-NH3 ice mix
    "HM80_rock": Di_mat_type["HM80"] * type_factor + 2,  # SiO2-MgO-FeS-FeO rock mix
    # SESAME
    "SESAME_iron": Di_mat_type["SESAME"] * type_factor,  # 2140
    "SESAME_basalt": Di_mat_type["SESAME"] * type_factor + 1,  # 7530
    "SESAME_water": Di_mat_type["SESAME"] * type_factor + 2,  # 7154
    "SS08_water": Di_mat_type["SESAME"] * type_factor + 3,  # Senft & Stewart (2008)
    "AQUA": Di_mat_type["SESAME"] * type_factor + 4,  # Haldemann+2020
    "CMS19_H": Di_mat_type["SESAME"] * type_factor + 5,  # Chabrier+2019 Hydrogen
    "CMS19_He": Di_mat_type["SESAME"] * type_factor + 6,  # Helium
    "CMS19_HHe": Di_mat_type["SESAME"] * type_factor + 7,  # H/He mixture Y=0.275
    # ANEOS
    "ANEOS_forsterite": Di_mat_type["ANEOS"] * type_factor,  # Stewart et al. (2019)
    "ANEOS_iron": Di_mat_type["ANEOS"] * type_factor + 1,  # Stewart (2020)
    "ANEOS_Fe85Si15": Di_mat_type["ANEOS"] * type_factor + 2,  # Stewart (2020)
}
# Invert so the ID are the keys
Di_id_mat = {mat_id: mat for mat, mat_id in Di_mat_id.items()}