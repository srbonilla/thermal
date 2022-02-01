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

# ID offset to distinguish different bodies e.g. impactor from target
id_body = 10000
A1_mat = [mat for mat in Di_mat_id.keys()]
A1_id = [id for id in Di_mat_id.values()]
for mat, id in zip(A1_mat, A1_id):
    Di_mat_id[mat + "_2"] = id + id_body

# Colours
Di_mat_colour = {
    # Ideal Gas
    "idg_HHe": "#1199ff",
    "idg_HHe_2": "#66ddee",
    "idg_N2": "#1199ff",
    "idg_N2_2": "#66ddee",
    "idg_CO2": "#1199ff",
    "idg_CO2_2": "#66ddee",
    # Tillotson
    "Til_iron": "#808080",
    "Til_iron_2": "#775533",
    "Til_granite": "#dd4400",
    "Til_granite_2": "#ffdd00",
    "Til_water": "#4169E1",
    "Til_water_2": "#4169E1",
    "Til_basalt": "#dd4400",
    "Til_basalt_2": "#ffdd00",
    # Hubbard & MacFarlane (1980) Uranus/Neptune
    "HM80_HHe": "#1199ff",
    "HM80_HHe_2": "#66ddee",
    "HM80_ice": "#B0C4DE",
    "HM80_ice_2": "#A080D0",
    "HM80_rock": "#708090",
    "HM80_rock_2": "#706050",
    # SESAME
    "SESAME_iron": "#808080",
    "SESAME_iron_2": "#775533",
    "SESAME_basalt": "#dd4400",
    "SESAME_basalt_2": "#ffdd00",
    "SESAME_water": "#4169E1",
    "SESAME_water_2": "#4169E1",
    "SS08_water": "#4169E1",
    "SS08_water_2": "#4169E1",
    "AQUA": "#4169E1",
    "AQUA_2": "#4169E1",
    "CMS19_H": "#1199ff",
    "CMS19_H_2": "#66ddee",
    "CMS19_He": "#1199ff",
    "CMS19_He_2": "#66ddee",
    "CMS19_HHe": "#1199ff",
    "CMS19_HHe_2": "#66ddee",
    # ANEOS
    "ANEOS_forsterite": "#dd4400",
    "ANEOS_forsterite_2": "#ffdd00",
    "ANEOS_iron": "#808080",
    "ANEOS_iron_2": "#775533",
    "ANEOS_Fe85Si15": "#808080",
    "ANEOS_Fe85Si15_2": "#775533",
}

# Invert so the ID are the keys
Di_id_type = {type_id: mat for mat, type_id in Di_mat_type.items()}
Di_id_colour = {Di_mat_id[mat]: c for mat, c in Di_mat_colour.items()}