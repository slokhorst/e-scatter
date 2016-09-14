from cstool.parse_input import read_input, pprint_settings, cstool_model


if __name__ == "__main__":
    s = read_input("silicondioxide.json")

    print(pprint_settings(cstool_model, s))
    print()
    print("Phonon loss: {:~P}".format(s.phonon_loss))
    print("Total molar weight: {:~P}".format(s.M_tot))
    print("Number density: {:~P}".format(s.rho_n))
