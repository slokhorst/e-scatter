from cstool.parse_input import cstool_model, read_input, pprint_settings

if __name__ == "__main__":
    s = pprint_settings(cstool_model, read_input("silicon.json"))
    print(s)

