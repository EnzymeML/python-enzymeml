import enzymeml.enzymeml as enzml
import enzymeml.tables as tables
import libsbml as sbml
import xlrd

#################################################################################
# This script works under following conditions at the moment:                   #
# The reaction is build in the following way:                                   #
# Substrate_1 + Substrate_2 (+Enzyme) -> Product_1 + Product_2 (+ Enzyme)       #
# <- and <-> are also handled.                                                  #
#                                                                               #
# Multiple reactions or enzymes are not beeing handled  and must be implemented #
# into the template file and also into this script.                             #
# Additionally the initial conditions are only entered if the exist at t=0.     #
#################################################################################

_unit_manager = dict()


# From BioCatNet_to_EnzymeML.py
def get_unit(enz, name):
    global _unit_manager

    if enz.id not in _unit_manager:
        _unit_manager[enz.id] = dict()

    units = _unit_manager[enz.id]

    if name not in units:
        if name == "%":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "percent",
                                               "units": [{"kind": sbml.UNIT_KIND_DIMENSIONLESS, "scale": -2}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000187", un)

            units["%"] = un
        elif name == "kelvin":
            units["kelvin"] = sbml.UNIT_KIND_KELVIN
        elif name == "ms":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "millisecond",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "scale": -3}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000028", un)

            units["ms"] = un
        elif name == "s":
            units["s"] = sbml.UNIT_KIND_SECOND
        elif name == "min":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "minute",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "multiplier": 60}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000031", un)

            units["min"] = un
        elif name == "h":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "hour",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "multiplier": 3600}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000032", un)

            units["h"] = un
        elif name == "nmol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "nmol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -9},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000062", un)

            units["nmol/l"] = un
        elif name == "umol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "umol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -6},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000064", un)

            units["umol/l"] = un
        elif name == "mmol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mmol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -3},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000065", un)

            units["mmol/l"] = un
        elif name == "mol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000062", un)

            units["mol/l"] = un
        elif name == "nmol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "nmol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -9}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000041", un)

            units["nmol"] = un
        elif name == "umol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "umol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -6}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000039", un)

            units["umol"] = un
        elif name == "mmol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mmol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -3}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000040", un)

            units["mmol"] = un
        elif name == "mol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000013", un)

            units["mol"] = un
        elif name == "rpm":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "revolution/min",
                                               "units": [{"kind": sbml.UNIT_KIND_DIMENSIONLESS},
                                                         {"kind": sbml.UNIT_KIND_SECOND, "exponent": -1,
                                                          "multiplier:": 60}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/NCIT:C70469", un)

            units["rpm"] = un
        elif name == "l":
            units["l"] = sbml.UNIT_KIND_LITRE
        elif name == "ml":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "ml",
                                               "units": [{"kind": sbml.UNIT_KIND_LITRE, "scale": -3}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000098", un)

            units["ml"] = un
        elif name == "ul":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "ul",
                                               "units": [{"kind": sbml.UNIT_KIND_LITRE, "scale": -6}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000101", un)

            units["ul"] = un
        elif name == "nl":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "nl",
                                               "units": [{"kind": sbml.UNIT_KIND_LITRE, "scale": -9}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000102", un)

            units["nl"] = un
        else:
            print("Unknown unit %s found. Illegal id is used." % name)
            units[name] = name

    return units[name]


def get_ranges(sheet, nameList):
    ranges = dict()

    name = None
    for row in range(0, sheet.nrows):
        cell = str(sheet.cell_value(row, 0))

        if cell.lower() in nameList:
            if name is not None:
                ranges[name].append(row - 1)
            name = cell.lower()
            ranges[name] = [row + 2]

    if name is not None:
        ranges[name].append(sheet.nrows - 1)

    return ranges


class Settings:
    def __init__(self):
        self.org = "Universität Stuttgart"


class UserInformationSheet:
    def __init__(self, sheet):
        self.user_name = sheet.cell_value(2, 0)
        self.user_email = sheet.cell_value(2, 1)

        self.exp_name = sheet.cell_value(6, 0)
        self.exp_set_name = sheet.cell_value(6, 1)


class ConditionsSheet:
    def __init__(self, sheet):
        self.cond_descr = sheet.cell_value(2, 0)
        self.cond_temp = sheet.cell_value(2, 1)
        self.cond_pH = sheet.cell_value(2, 2)
        self.cond_pressure = sheet.cell_value(4, 0)
        self.cond_pressure_unit = sheet.cell_value(4, 1)
        self.cond_shaking = sheet.cell_value(4, 2)
        self.cond_init_volume = sheet.cell_value(6, 0)
        self.cond_init_volume_unit = sheet.cell_value(6, 1)

        # Only one enzyme/sequence
        self.enzyme_name = sheet.cell_value(10, 0)
        self.enzyme_expr_host = sheet.cell_value(10, 1)
        self.enzyme_prep = sheet.cell_value(10, 2)
        self.enzyme_application = sheet.cell_value(10, 3)

        # Buffer
        self.buffer_name = sheet.cell_value(14, 0)
        self.buffer_descr = sheet.cell_value(14, 1)
        self.buffer_conc = sheet.cell_value(16, 0)
        self.buffer_conc_unit = sheet.cell_value(16, 1)
        self.buffer_additives = sheet.cell_value(16, 2)


class ReactionSheet:
    class Compound:
        def __init__(self, name, smiles, prep):
            self.name = name
            self.smiles = smiles
            self.preparation = prep

        def has_smiles(self):
            return self.smiles is not None

        def has_notes(self):
            return self.preparation is not None

    class Reaction:
        def __init__(self, name, reac):
            self.name = name
            self.reaction = reac

    class Species:
        def __init__(self, measure, name, tv, tu, cv, cu, av, au, vv, vu):
            self.measurement_no = measure
            self.name = name
            self.time_value = tv
            self.time_unit = tu

            self.conc_val = cv
            self.conc_unit = cu
            self.amount_val = av
            self.amount_unit = au
            self.volume_val = vv
            self.volume_unit = vu

        def has_conc(self):
            return self.conc_val is not None and self.conc_unit is not None

        def has_amount(self):
            return self.amount_val is not None and self.amount_unit is not None

        def has_volume(self):
            return self.volume_val is not None and self.volume_unit is not None

    def __init__(self, sheet):  # Reactions should ALWAYS be 1, but infrastructure for multiple enzymes is there
        self.compounds = dict()
        self.reactions = dict()
        self.enzymes = dict()
        self.substrates = dict()
        self.additives = dict()

        ranges = get_ranges(sheet, ["compounds", "reaction", "enzymes", "substrates", "additives", "comments"])

        r = ranges["compounds"]
        for i in range(r[0], r[1]):
            c0 = sheet.cell_value(i, 0)
            if c0 != "":
                c = ReactionSheet.Compound(
                    c0, sheet.cell_value(i, 1), sheet.cell_value(i, 2)
                )
                self.compounds[c.name] = c

        r = ranges["reaction"]
        for i in range(r[0], r[1]):
            c0 = sheet.cell_value(i, 0)
            if c0 != "":
                re = ReactionSheet.Reaction(c0, sheet.cell_value(i, 1))
                self.reactions[re.name] = re

        r = ranges["enzymes"]
        for i in range(r[0], r[1]):
            c1 = sheet.cell_value(i, 1)
            if c1 != "":
                sp = ReactionSheet.Species(
                    sheet.cell_value(i, 0),
                    c1,
                    sheet.cell_value(i, 2),
                    sheet.cell_value(i, 3),
                    sheet.cell_value(i, 4),
                    sheet.cell_value(i, 5),
                    sheet.cell_value(i, 6),
                    sheet.cell_value(i, 7),
                    sheet.cell_value(i, 8),
                    sheet.cell_value(i, 9)
                )
                if sp.name not in self.enzymes:
                    self.enzymes[sp.name] = list()
                self.enzymes[sp.name].append(sp)

        r = ranges["substrates"]
        for i in range(r[0], r[1]):
            c1 = sheet.cell_value(i, 1)
            if c1 != "":
                sp = ReactionSheet.Species(
                    sheet.cell_value(i, 0),
                    c1,
                    sheet.cell_value(i, 2),
                    sheet.cell_value(i, 3),
                    sheet.cell_value(i, 4),
                    sheet.cell_value(i, 5),
                    sheet.cell_value(i, 6),
                    sheet.cell_value(i, 7),
                    sheet.cell_value(i, 8),
                    sheet.cell_value(i, 9)
                )
                if sp.name not in self.substrates:
                    self.substrates[sp.name] = list()
                self.substrates[sp.name].append(sp)

        r = ranges["additives"]
        for i in range(r[0], r[1]):
            c1 = sheet.cell_value(i, 1)
            if c1 != "":
                sp = ReactionSheet.Species(
                    sheet.cell_value(i, 0),
                    c1,
                    sheet.cell_value(i, 2),
                    sheet.cell_value(i, 3),
                    sheet.cell_value(i, 4),
                    sheet.cell_value(i, 5),
                    sheet.cell_value(i, 6),
                    sheet.cell_value(i, 7),
                    sheet.cell_value(i, 8),
                    sheet.cell_value(i, 9)
                )
                if sp.name not in self.additives:
                    self.additives[sp.name] = list()
                self.additives[sp.name].append(sp)

    def get_substrates(self, reac_name):
        reac = self.reactions[reac_name]

        if reac is None:
            return None

        substrate_str = None

        if "->" in reac.reaction:
            substrate_str = reac.reaction.split("->")[0]
        elif "<-" in reac.reaction:
            substrate_str = reac.reaction.split("<-")[1]
        elif "<->" in reac.reaction:
            substrate_str = reac.reaction.split("<->")[0]
        else:
            raise RuntimeError("Unidentified reaction found. Could not find the equation arrows" +
                               "'->', '<-' or '<->' in '%s'" % reac.reaction)

        substrates = substrate_str.split("+")
        subs_list = list()

        for s in substrates:
            s_split = s.split()  # Removes spaces

            stochio = 1

            if len(s_split) > 2:
                raise RuntimeError("Reaction must be displayed with only a number and the name.")
            elif len(s_split) == 2:
                s = s_split[1]
                stochio = int(s_split[0])
            elif len(s_split) == 1:
                s = s_split[0]
            else:
                raise RuntimeError("Missing component in the reaction.")

            if not self.is_enzyme(s):
                subs_list.append((stochio, s))

        return subs_list

    def get_products(self, reac_name):
        reac = self.reactions[reac_name]

        if reac is None:
            return None

        product_str = None

        if "->" in reac.reaction:
            product_str = reac.reaction.split("->")[1]
        elif "<-" in reac.reaction:
            product_str = reac.reaction.split("<-")[0]
        elif "<->" in reac.reaction:
            product_str = reac.reaction.split("<->")[1]
        else:
            raise RuntimeError("Unidentified reaction found. Could not find the equation arrows" +
                               "'->', '<-' or '<->' in '%s'" % reac.reaction)

        products = product_str.split("+")
        prods_list = list()

        for p in products:
            p_split = p.split()  # Removes spaces
            stochio = 1

            if len(p_split) > 2:
                raise RuntimeError("Reaction must be displayed with only a number and the name.")
            elif len(p_split) == 2:
                p = p_split[1]
                stochio = int(p_split[0])
            elif len(p_split) == 1:
                p = p_split[0]
            else:
                raise RuntimeError("Missing component in the reaction.")

            if not self.is_enzyme(p):
                prods_list.append((stochio, p))

        return prods_list

    def is_enzyme(self, name):
        if type(name) is not str:
            raise TypeError("The name is not a string.")
        for e in self.enzymes:
            if e == name:
                return True

        return False

    # Return None, if not found, else:
    # (is_conc, amount, unit)
    def get_init_species(self, name):
        sl = None
        if name in self.substrates:
            sl = self.substrates[name]
        elif name in self.enzymes:
            sl = self.enzymes[name]
        elif name in self.additives:
            sl = self.additives[name]
        else:
            return None

        for s in sl:
            if s.time_value == 0:
                continue

            if s.has_conc():
                return True, s.conc_val, s.conc_unit
            elif s.has_amount():
                return False, s.amount_val, s.amount_unit
            elif s.has_volume():
                return False, s.volume_val, s.volume_unit

        return None


class SequencesSheet:
    class Sequence:
        def __init__(self, name, src, seq):
            self.name = name
            self.src_org = src
            self.sequence = seq

    def __init__(self, sheet):
        self.sequences = dict()

        ranges = get_ranges(sheet, ["sequences", "comments"])

        r = ranges["sequences"]
        for i in range(r[0], r[1]):
            c0 = sheet.cell_value(i, 0)
            if c0 != "":
                s = SequencesSheet.Sequence(
                    c0,
                    sheet.cell_value(i, 1),
                    sheet.cell_value(i, 2)
                )
                self.sequences[c0] = s


class MeasurementSheet:
    class Measurement:
        def __init__(self, mno, rno, name, tv, tu, cv, cu, method):
            self.measurement_num = mno
            self.replication_num = rno
            self.name = name
            self.time_val = tv
            self.time_unit = tu
            self.conc_val = cv
            self.conc_unit = cu
            self.method = method

    def __init__(self, sheet):
        self.measurements = dict()
        self.measurement_numbers = set()

        ranges = get_ranges(sheet, ["measured compounds", "comments"])

        r = ranges["measured compounds"]
        for i in range(r[0], r[1]):
            c2 = sheet.cell_value(i, 2)
            if c2 != "":
                if c2 not in self.measurements:
                    self.measurements[c2] = list()
                mn = sheet.cell_value(i, 0)
                m = MeasurementSheet.Measurement(
                    mn,
                    sheet.cell_value(i, 1),
                    c2,
                    sheet.cell_value(i, 3),
                    sheet.cell_value(i, 4),
                    sheet.cell_value(i, 5),
                    sheet.cell_value(i, 6),
                    sheet.cell_value(i, 7)  # This one is not used at the moment
                )
                if c2 not in self.measurements:
                    self.measurements[c2] = list()
                self.measurements[c2].append(m)
                self.measurement_numbers.add(mn)

    def get_replications(self, name, measure_no):
        measures = self.measurements[name]

        out = dict()

        for measure in measures:
            if measure.measurement_num != measure_no:
                continue
            if measure.replication_num not in out:
                out[measure.replication_num] = list()
            out[measure.replication_num].append(measure)

        return out


def parse(name, enz, sett):
    wb = xlrd.open_workbook(name)

    #########################
    # Load user information #
    #########################
    uis = wb.sheet_by_name("user information")
    tuis = UserInformationSheet(uis)

    enz.add_creator(tuis.user_name, tuis.user_name, tuis.user_email, sett.org)
    enz.add(enzml.key.MAIN_META_CREATOR, {
        "family": tuis.user_name, "given": tuis.user_name, "email": tuis.user_email, "org": sett.org}
            )

    ###################
    # Load conditions #
    ###################
    cs = wb.sheet_by_name("conditions")
    tcs = ConditionsSheet(cs)

    comp = enz.add(enzml.key.MAIN_COMPARTMENT, {
        "size": tcs.cond_init_volume, "units": get_unit(enz, tcs.cond_init_volume_unit), "constant": True
    })

    ########################
    # Initialize Reactions #
    ########################
    rs = wb.sheet_by_name("reaction")
    trs = ReactionSheet(rs)

    reac_ids = dict()
    for r in trs.reactions:
        re = trs.reactions[r]
        reac_ids[r] = enz.add(enzml.key.MAIN_REACTION, {"name": r, "reversible": "<->" in re.reaction})

    ##################
    # Load sequences #
    ##################

    ss = wb.sheet_by_name("sequences")
    tss = SequencesSheet(ss)

    ################
    # Load species #
    ################
    species = dict()
    subs = list()
    prods = list()
    for reac in reac_ids:
        subs += trs.get_substrates(reac)
        prods += trs.get_products(reac)

    for s in trs.enzymes:
        # sp = trs.enzymes[s]
        obj = {"name": s, "compartment": comp, "type": enzml.ontology.SBO_METABOLITE, "constant": True}

        init = trs.get_init_species(s)
        if init is not None:
            if init[0]:
                obj["init_conc"] = init[1]
            else:
                obj["init_amount"] = init[1]

            obj["units"] = get_unit(enz, init[2])

        sp_id = enz.add(enzml.key.MAIN_SPECIES, obj)

        c = trs.enzymes[s]

        notes = ""
        if tcs.enzyme_expr_host != "":
            notes += "Expression host: %s" % tcs.enzyme_expr_host

        if tcs.enzyme_prep != "":
            if notes != "":
                notes += "<br/>"
            notes += "Preparation: %s" % tcs.enzyme_prep

        if tcs.enzyme_application != "":
            if notes != "":
                notes += "<br/>"
            notes += "Application: %s" % tcs.enzyme_application

        if notes != "":
            enz.add(enzml.key.UNSPECIFIC_NOTE, notes, sp_id)

        seqs = tss.sequences[s]
        enz.add(enzml.key.MAIN_SPECIES_PROTEIN,
                {
                    "sequence": seqs.sequence,
                    "hasTaxon": "taxonomy:%s" % seqs.src_org, "occursIn": "taxonomy:%s" % tcs.enzyme_expr_host
                }, sp_id)

        species[s] = sp_id

    metabolits = set(subs) & set(prods)
    for s in metabolits:
        obj = {"name": s, "compartment": comp, "type": enzml.ontology.SBO_METABOLITE, "constant": False}

        init = s.get_init_species(s)
        if init is not None:
            if init[0]:
                obj["init_conc"] = init[1]
            else:
                obj["init_amount"] = init[1]

            obj["units"] = get_unit(enz, init[2])

        sp_id = enz.add(enzml.key.MAIN_SPECIES, obj)

        c = trs.compounds[s]
        if c.has_smiles():
            enz.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": c.smiles})
        if c.has_notes():
            enz.add(enzml.key.UNSPECIFIC_NOTE, c.preparation, sp_id)

        species[s] = sp_id

    for s in subs:
        if s not in metabolits:
            obj = {
                "name": s[1], "compartment": comp, "type": enzml.ontology.SBO_SUBSTRATE,
                "constant": False
            }

            init = trs.get_init_species(s)
            if init is not None:
                if init[0]:
                    obj["init_conc"] = init[1]
                else:
                    obj["init_amount"] = init[1]

                obj["units"] = get_unit(enz, init[2])

            sp_id = enz.add(enzml.key.MAIN_SPECIES, obj)

            c = trs.compounds[s[1]]
            if c.has_smiles():
                enz.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": c.smiles}, sp_id)
            if c.has_notes():
                enz.add(enzml.key.UNSPECIFIC_NOTE, c.preparation, sp_id)

            species[s[1]] = sp_id

    for s in prods:
        if s not in metabolits:
            obj = {
                "name": s[1], "compartment": comp, "type": enzml.ontology.SBO_PRODUCT,
                "constant": False
            }

            init = trs.get_init_species(s)
            if init is not None:
                if init[0]:
                    obj["init_conc"] = init[1]
                else:
                    obj["init_amount"] = init[1]

                obj["units"] = get_unit(enz, init[2])

            sp_id = enz.add(enzml.key.MAIN_SPECIES, obj)

            c = trs.compounds[s[1]]
            if c.has_smiles():
                enz.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": c.smiles}, sp_id)
            if c.has_notes():
                enz.add(enzml.key.UNSPECIFIC_NOTE, c.preparation, sp_id)

            species[s[1]] = sp_id

    ###################
    # Setup reactions #
    ###################
    for reac in reac_ids:
        subs = trs.get_substrates(reac)
        for s in subs:
            enz.add(enzml.key.MODEL_REACTION_REACTANTS, {"id": species[s[1]], "stochiometry": s[0]}, reac_ids[reac])

        prods = trs.get_products(reac)
        for p in prods:
            enz.add(enzml.key.MODEL_REACTION_PRODUCTS, {"id": species[p[1]], "stochiometry": p[0]}, reac_ids[reac])

        conds = dict()
        if tcs.cond_temp != "":
            conds["temperature"] = (float(tcs.cond_temp) + 273.15, get_unit(enz, "kelvin"))

        if tcs.cond_pressure != "":
            conds["pressure"] = (float(tcs.cond_pressure), get_unit(enz, tcs.cond_pressure_unit))

        if tcs.cond_pH != "":
            conds["ph"] = float(tcs.cond_pH)

        if tcs.cond_shaking != "":
            conds["shaking"] = (float(tcs.cond_shaking), get_unit(enz, "rpm"))

        enz.add(enzml.key.MAIN_REACTION_CONDITION, conds, reac_ids[reac])

    ######################
    # Writing CSV files. #
    ######################
    mcs = wb.sheet_by_name("measured compounds")
    tmcs = MeasurementSheet(mcs)

    f = enzml.EnzymeMLFormat()
    enz.add(enzml.key.MAIN_DATA_FORMAT, f)
    tc = enzml.create_column(enzml.COLUMN_TYPE_TIME, "s")
    f.add_column(tc)

    cdt = tables.ColumnDependentTable("time")
    for compound in tmcs.measurements:
        replica = tmcs.get_replications(compound, 1)
        for repl_name in replica:
            r = replica[repl_name]
            col_name = "%s_%s" % (compound, repl_name)
            cdt.init_column(col_name)
            c = enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, species[r[0].name], get_unit(enz, r[0].conc_unit))
            f.add_column(c)
            for m in r:
                cdt.add(m.time_val, col_name, m.conc_val)

    csv = enzml.EnzymeMLCSV(f, name="data")
    enz.add_csv(csv)
    csv_file = enz.add(enzml.key.MAIN_DATA_FILE, {"file": csv.location, "format": f.sid})

    measure = enz.add(enzml.key.MAIN_DATA_MEASUREMENTS,
                      {
                          "file": csv_file, "start": 1, "stop": cdt.nrows(),
                          "name": "Experiment measurement data."
                      })

    known_repl = set()
    for i in range(1, len(f.columns)):
        col = f.columns[i]
        if col.replica in known_repl:
            continue
        known_repl.add(col.replica)
        repl_enzml = enzml.EnzymeMLReplica(measure, col.replica)
        enz.add(enzml.key.MAIN_REACTION_REPLICAS, repl_enzml, list(reac_ids.values())[0])

    columns = cdt.as_columns()

    csv.add_column(columns[0])

    for i in range(1, len(columns)):
        csv.add_column(columns[i])

    pass


if __name__ == "__main__":
    print("Parsing templates/test_template.xlsx File to EnzymeML")
    enz = enzml.EnzymeML("ParsedTemplate")
    settings = Settings()
    parse("templates/test_template.xlsx", enz, settings)

    print("Checking for completeness...")
    # Do something here

    print("Creating archive...")
    enz.create_archive()
    print("Finished parsing.")
