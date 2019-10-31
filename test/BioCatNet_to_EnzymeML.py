# --------------------------------------------------------------------------------------------------
# Read from/ write to Firebird Database via fdb
# On basis of Jens Lehmans script
# --------------------------------------------------------------------------------------------------

import libsbml as sbml
import fdb
import enzymeml.enzymeml as enzml
import statistics
import enzymeml.tables as tables
import traceback


# Used to handle replications for writing the CSV files
class Replication:
    class Entry:
        def __init__(self, time, conc):
            self.time = time
            self.conc = conc

        def __lt__(self, other):
            return self.time < other.time

        def __str__(self):
            return "Time %s: Concentration %s" % (self.time, self.conc)

    def __init__(self, data):
        self.entries = list()
        self.time_unit = None
        self.conc_unit = None

        for dat in data:
            if dat[1] != self.time_unit:
                if self.time_unit is None:
                    self.time_unit = dat[1]
                else:
                    raise RuntimeError("Inconsistent Time Units data with Time %s and unit %s" % (dat[0], dat[1]))
            if dat[3] != self.conc_unit:
                if self.conc_unit is None:
                    self.conc_unit = dat[3]
                else:
                    raise RuntimeError(
                        "Inconsistent Conc Units data with Time/Conc %s/%s and unit %s <-> %s"
                        % (dat[0], dat[2], dat[3], self.conc_unit))

            self.entries.append(Replication.Entry(dat[0], dat[2]))

        self.entries.sort()

    def initial_conc(self):
        return self.entries[0].conc

    def add_to_table(self, table, name):
        table.init_column(name)

        for e in self.entries:
            table.add(e.time, name, e.conc)


class ReplicationSet:
    def __init__(self):
        self.replica = list()

    def add(self, repl):
        self.replica.append(repl)

    def get_time_unti(self):
        return self.replica[-1].time_unit

    def get_conc_unit(self):
        return self.replica[-1].conc_unit

    def print_validation(self):
        le = len(self.replica)
        if le < 2:
            return True

        tu = self.replica[0].time_unit
        cu = self.replica[0].conc_unit

        time_stamps = list()
        for e in self.replica[0].entries:
            time_stamps.append(e.time)

        if len(time_stamps) == 0:
            time_stamps.append(0)  # To prevent errors if nothing is in the first one for whatever reason

        for i in range(1, le):
            repl = self.replica[i]

            if tu != repl.time_unit or cu != repl.conc_unit:
                print("The replications contain different units and cannot be saved.")
                return False

            for e in repl.entries:
                if e.time not in time_stamps:
                    if e.time > time_stamps[-1]:
                        time_stamps.append(e.time)
                    else:
                        print("The replication %i has different time stamps than the ones before." % i)
                        return False

        return True

    def __len__(self):
        return self.replica.__len__()


def init_conc_mean(replica):
    inits = list()
    for rep in replica.replica:
        inits.append(rep.initial_conc())

    mean = statistics.mean(inits)

    stdev = 0
    if len(inits) > 1:
        stdev = statistics.stdev(inits)

    return mean, stdev


_unit_manager = dict()


def _get_unit_descr(name):
    descr_str = """
    SELECT DESCRIPTION FROM UNITS WHERE NAME = '%s'
    """ % name

    return cur.execute(descr_str).fetchone()[-1]


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
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("%"), un)

            units["%"] = un
        elif name == "kelvin":
            units["kelvin"] = sbml.UNIT_KIND_KELVIN
        elif name == "ms":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "millisecond",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "scale": -3}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000028", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("ms"), un)

            units["ms"] = un
        elif name == "s":
            units["s"] = sbml.UNIT_KIND_SECOND
        elif name == "min":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "minute",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "multiplier": 60}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000031", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("min"), un)

            units["min"] = un
        elif name == "h":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "hour",
                                               "units": [{"kind": sbml.UNIT_KIND_SECOND, "multiplier": 3600}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000032", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("h"), un)

            units["h"] = un
        elif name == "nmol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "nmol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -9},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000062", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("nmol/l"), un)

            units["nmol/l"] = un
        elif name == "umol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "umol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -6},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000064", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("umol/l"), un)

            units["umol/l"] = un
        elif name == "mmol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mmol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -3},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000065", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("mmol/l"), un)

            units["mmol/l"] = un
        elif name == "mol/l":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mol/l",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE},
                                                         {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000062", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("mol/l"), un)

            units["mol/l"] = un
        elif name == "nmol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "nmol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -9}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000041", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("nmol"), un)

            units["nmol"] = un
        elif name == "umol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "umol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -6}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000039", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("umol"), un)

            units["umol"] = un
        elif name == "mmol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mmol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -3}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000040", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("mmol"), un)

            units["mmol"] = un
        elif name == "mol":
            un = enz.add(enzml.key.MAIN_UNIT, {"name": "mol",
                                               "units": [{"kind": sbml.UNIT_KIND_MOLE}]
                                               })
            enz.add(enzml.key.MAIN_UNIT_IS, "https://identifiers.org/UO:0000013", un)
            enz.add(enzml.key.UNSPECIFIC_NOTE, _get_unit_descr("mol"), un)

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


def load_model(cur, exp_id):
    reac_id_str = """
    SELECT REACTION_ID FROM EXP_RKT_LINK
    WHERE EXPERIMENT_ID = %s;""" % exp_id
    # JOIN EXPERIMENTS ON EXP_RKT_LINK.EXPERIMENT_ID = EXPERIMENTS.EXPERIMENT_ID
    reac_id = cur.execute(reac_id_str).fetchone()[-1]

    cond_str = """
    SELECT USER_ID, CONDITION_ID, DESCRIPTION, SUBMIT_DATE, GROUP_ID FROM EXPERIMENTS
    WHERE EXPERIMENT_ID = %s;""" % exp_id
    cond_id = cur.execute(cond_str).fetchall()[-1]

    enz = enzml.EnzymeML("BioCatNet-Test")

    # Species definition
    s_subs = """
                SELECT NAME, REACTION_COMPOUNDS.COMPOUND_ID FROM COMPOUND_NAMES
                JOIN REACTION_COMPOUNDS ON COMPOUND_NAMES.COMPOUND_ID = REACTION_COMPOUNDS.COMPOUND_ID 
                WHERE REACTION_COMPOUNDS.REACTION_ID = %s AND PRODUCT = 0;""" % reac_id
    substrates = cur.execute(s_subs).fetchall()

    s_prod = """
                SELECT NAME, REACTION_COMPOUNDS.COMPOUND_ID FROM COMPOUND_NAMES
                JOIN REACTION_COMPOUNDS ON COMPOUND_NAMES.COMPOUND_ID = REACTION_COMPOUNDS.COMPOUND_ID 
                WHERE REACTION_COMPOUNDS.REACTION_ID = %s AND PRODUCT = 1;""" % reac_id
    products = cur.execute(s_prod).fetchall()

    # Searching for initial Concentrations
    data_str = """
                SELECT CONCENTRATION FROM COMPOUND_MEASUREMENTS
                WHERE EXPERIMENT_ID = %s AND POINTINTIME = 0;""" % exp_id
    data = cur.execute(data_str).fetchall()

    replic_str = """SELECT DISTINCT COMPOUND_MEASUREMENTS.REPLICATION
            FROM COMPOUND_MEASUREMENTS WHERE EXPERIMENT_ID = %s""" % exp_id
    replica = cur.execute(replic_str).fetchall()

    reac_desc_str = """SELECT DESCRIPTION FROM REACTIONS WHERE REACTION_ID = %s""" % reac_id
    reac_desc = cur.execute(reac_desc_str).fetchone()[-1]

    # Create compartment
    comp = enz.add(enzml.key.MAIN_COMPARTMENT, {"dimensions": 3, "constant": True})

    # Create reaction
    reac = enz.add(enzml.key.MAIN_REACTION, {"name": reac_desc, "reversible": False})

    unit_ids = dict()

    def unit_from_id(unit):
        if unit not in unit_ids:
            u_str = """
            SELECT NAME FROM UNITS WHERE UNIT_ID = %s
            """ % unit
            unit_ids[unit] = cur.execute(u_str).fetchone()[-1]
        return unit_ids[unit]

    form = enzml.EnzymeMLFormat()
    form.add_column(enzml.create_column(enzml.COLUMN_TYPE_TIME, "seconds"))

    enz.add(enzml.key.MAIN_DATA_FORMAT, form)

    csv = enzml.EnzymeMLCSV(form, name="data")
    enz.add_csv(csv)
    csv_file = enz.add(enzml.key.MAIN_DATA_FILE, {"file": csv.location, "format": form.sid})

    deptable = tables.ColumnDependentTable("time")

    # Adding substrates
    for subs in substrates:
        # subs: (Name, ID)
        subs_name = subs[0]

        replications = ReplicationSet()

        for repl in replica:
            repl_data_str = """
            SELECT COMPOUND_MEASUREMENTS.POINTINTIME, COMPOUND_MEASUREMENTS.TIME_UNIT,
                COMPOUND_MEASUREMENTS.CONCENTRATION, COMPOUND_MEASUREMENTS.CONC_UNIT
            FROM COMPOUND_MEASUREMENTS
            WHERE EXPERIMENT_ID = %s AND COMPOUND_ID = %s AND REPLICATION = %s""" % (exp_id, subs[1], repl[-1])
            data = cur.execute(repl_data_str).fetchall()
            # data: Time, Time-Unit-ID, Concentration, Conc-Unit-ID

            replications.add(Replication(data))

        mean, stdev = init_conc_mean(replications)

        obj = {"name": subs_name, "compartment": comp,
               "type": enzml.ontology.SBO_SUBSTRATE, "constant": False,
               "init_conc": mean, "units": get_unit(enz, unit_from_id(replications.get_conc_unit()))
               }

        if stdev > 0:
            obj["stdev"] = stdev

        sub = enz.add(enzml.key.MAIN_SPECIES, obj)
        enz.add(enzml.key.MAIN_REACTION_REACTANTS, {"id": sub}, reac)

        smiles_data_str = """
        SELECT SMILES FROM COMPOUNDS WHERE COMPOUND_ID = %s
        """ % subs[1]
        smiles = cur.execute(smiles_data_str).fetchone()
        if len(smiles) > 0 and smiles[-1] is not None:
            enz.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": smiles[-1]}, sub)
        else:
            print("Could not find smiles code for %s with the name '%s'." % (subs[1], subs_name))

        if replications.print_validation():
            for repl in replications.replica:
                col = enzml.create_column(
                    enzml.COLUMN_TYPE_CONCENTRATION, sub, get_unit(enz, unit_from_id(repl.conc_unit))
                )
                form.add_column(col)
                repl.add_to_table(deptable, "%s_%s" % (sub[0], col.replica))
                measure = enz.add(enzml.key.MAIN_DATA_MEASUREMENTS,
                                  {
                                      "file": csv_file, "start": 1, "stop": len(repl.entries),
                                      "name": "Measurement %s" % subs_name
                                  })
                repl_enzml = enzml.EnzymeMLReplica(measure, col.replica)
                enz.add(enzml.key.MAIN_REACTION_REPLICAS, repl_enzml, reac)

    # Adding products
    for prods in products:
        # prods: (Name, ID)
        prods_name = prods[0]

        replications = ReplicationSet()

        for repl in replica:
            repl_data_str = """
            SELECT COMPOUND_MEASUREMENTS.POINTINTIME, COMPOUND_MEASUREMENTS.TIME_UNIT,
                COMPOUND_MEASUREMENTS.CONCENTRATION, COMPOUND_MEASUREMENTS.CONC_UNIT
            FROM COMPOUND_MEASUREMENTS
            WHERE EXPERIMENT_ID = %s AND COMPOUND_ID = %s AND REPLICATION = %s""" % (exp_id, prods[1], repl[-1])
            data = cur.execute(repl_data_str).fetchall()
            # data: Time, Time-Unit-ID, Concentration, Conc-Unit-ID

            if len(data) > 0:
                replications.append(Replication(data))

        mean = 0
        stdev = -1
        units = get_unit(enz, "mol/l")

        if len(replications) > 0:
            mean, stdev = init_conc_mean(replications)
            units = get_unit(enz, unit_from_id(replications[-1].conc_unit))

        obj = {"name": prods_name, "compartment": comp,
               "type": enzml.ontology.SBO_PRODUCT, "constant": False,
               "init_conc": mean, "units": units}

        if not stdev < 0:
            obj["stdev"] = stdev

        prod = enz.add(enzml.key.MAIN_SPECIES, obj)

        enz.add(enzml.key.MAIN_REACTION_PRODUCTS, {"id": prod}, reac)

        smiles_data_str = """
                SELECT SMILES FROM COMPOUNDS WHERE COMPOUND_ID = %s
                """ % prods[1]
        smiles = cur.execute(smiles_data_str).fetchone()
        if len(smiles) > 0 and smiles[-1] is not None:
            enz.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": smiles[-1]}, prod)
        else:
            print("Could not find smiles code for %s with the name '%s'." % (prods[1], prods_name))

        if replications.print_validation():
            for repl in replications.replica:
                col = enzml.create_column(
                    enzml.COLUMN_TYPE_CONCENTRATION, sub, get_unit(enz, unit_from_id(repl.conc_unit))
                )
                form.add_column(col)
                repl.add_to_table(deptable, "%s_%s" % (sub[0], col.replica))
                measure = enz.add(enzml.key.MAIN_DATA_MEASUREMENTS,
                                  {
                                      "file": csv_file, "start": 1, "stop": len(repl.entries),
                                      "name": "Measurement %s" % subs_name
                                  })
                repl_enzml = enzml.EnzymeMLReplica(measure, col.replica)
                enz.add(enzml.key.MAIN_REACTION_REPLICAS, repl_enzml, reac)

    # Add all replication data saved in the time dependent table in the csv file
    for col in deptable.as_columns():
        csv.add_column(col)

    # Adding the enzyme
    enzyme_str = """
    SELECT PROTEINS.PROTEIN_NAME, ENZYME_FEEDS.*, SEQUENCES.SEQUENCE_ID FROM PROTEINS
    JOIN SEQUENCES ON PROTEINS.PROTEIN_ID = SEQUENCES.PROTEIN_ID
    JOIN ENZYME_FEEDS ON SEQUENCES.SEQUENCE_ID = ENZYME_FEEDS.SEQUENCE_ID
    WHERE ENZYME_FEEDS.EXPERIMENT_ID = %s;""" % exp_id

    # enzyme: name, enzyme_feed, sequence_id, units, preparation, expression_host, exp_id, time, amount, volume,
    #         time unit, amount unit, volume unit, conc unit, conc, seq_id
    enzyme = cur.execute(enzyme_str).fetchall()[-1]

    obj = {"name": enzyme[0], "compartment": comp,
           "type": enzml.ontology.SBO_ENZYME, "constant": True}

    if enzyme[14] is not None:  # Concentration
        obj["init_conc"] = enzyme[14]
        obj["units"] = get_unit(enz, unit_from_id(enzyme[13]))
    elif enzyme[8] is not None:  # Amount
        obj["init_conc"] = enzyme[8]
        obj["units"] = get_unit(enz, unit_from_id(enzyme[11]))
    elif enzyme[9] is not None:  # Volume
        obj["init_conc"] = enzyme[9]
        obj["units"] = get_unit(enz, unit_from_id(enzyme[12]))
    else:  # Nothing given
        obj["init_conc"] = 0
        obj["units"] = get_unit(enz, "mol/l")

    enz_id = enz.add(enzml.key.MAIN_SPECIES, obj)
    enz.add(enzml.key.MAIN_REACTION_MODIFIERS, {"id": enz_id}, reac)

    ec_str = """
    SELECT EC, PROTEINS.DESCRIPTION FROM PROTEINS
    JOIN SEQUENCES ON PROTEINS.PROTEIN_ID = SEQUENCES.PROTEIN_ID
    JOIN ENZYME_FEEDS ON SEQUENCES.SEQUENCE_ID = ENZYME_FEEDS.SEQUENCE_ID
    WHERE ENZYME_FEEDS.EXPERIMENT_ID = %s;""" % exp_id
    ec = cur.execute(ec_str).fetchall()[-1]

    enz.add(enzml.key.MAIN_REACTION_EC_CODE, ec[0], reac)
    enz.add(enzml.key.UNSPECIFIC_NOTE, ec[1], enz_id)

    seq_str = """
    SELECT LIST(AA,'') FROM POSITIONS
    WHERE SEQUENCE_ID = %s
    GROUP BY SEQUENCE_ID
    """ % enzyme[15]

    cur.execute(seq_str)
    sequence = cur.fetchone()[-1]

    enz.add(enzml.key.MAIN_SPECIES_PROTEIN, {"sequence": sequence}, enz_id)

    # Addition of Buffer
    cond_str = """
        SELECT * FROM CONDITIONS WHERE CONDITION_ID = %s
        """ % cond_id[1]
    # conds: id, buffer, description, shaking, temperature, ph, pressure, scale, duration, duration unit, pressure unit
    #        scale unit
    conds = cur.execute(cond_str).fetchall()[-1]

    buffer_str = """
    SELECT * FROM BUFFERS
    WHERE BUFFERS.BUFFER_ID = %s;""" % conds[1]
    # buffer: id, description, ions, conc, conc unit, name
    buffer = cur.execute(buffer_str).fetchall()[-1]

    sp_buffer = enz.add(enzml.key.MAIN_SPECIES, {"name": buffer[5], "compartment": comp,
                                                 "type": enzml.ontology.SBO_NEUTRAL_PARTICIPANT, "constant": True,
                                                 "init_conc": buffer[3],
                                                 "units": get_unit(enz, unit_from_id(buffer[4]))})
    enz.add(enzml.key.UNSPECIFIC_NOTE, buffer[1], sp_buffer)
    enz.add(enzml.key.MAIN_REACTION_MODIFIERS, {"id": sp_buffer}, reac)

    # Addition of conditions
    enz.add(enzml.key.MAIN_REACTION_CONDITION, {
        "ph": conds[5], "temperature": (float(conds[4]) + 273.13, sbml.UNIT_KIND_KELVIN),  # To kelvin
        "pressure": (conds[6], unit_from_id(conds[10])), "shaking": (conds[3], get_unit(enz, "rpm"))
    }, reac)

    # Add creator information
    enz.add_creator("Halupczok", "Colin", None, "Universität Stuttgart")  # FIXME Me as the creator of the archive

    user_str = "SELECT TITLE, FIRSTNAME, LASTNAME, EMAIL FROM USERS WHERE USER_ID = %s" % cond_id[0]
    user = cur.execute(user_str).fetchall()[-1]

    """
    For organisation
    SELECT GROUPS.GROUP_NAME FROM GROUPS
JOIN USER_GROUP ON GROUPS.GROUP_ID = USER_GROUP.GROUP_ID
WHERE USER_GROUP.USER_ID = 40
    """

    enz.add(enzml.key.MAIN_META_CREATOR, {"family": user[2], "given": user[1], "email": user[3]})
    date = cond_id[3]
    sbml_date = sbml.Date(date.year, date.month, date.day, date.hour, date.minute, date.second)
    enz.add(enzml.key.MAIN_META_DATES_CREATE, sbml_date)

    return enz


if __name__ == '__main__':
    # Building the connection and the cursor
    con = fdb.connect(dsn='URI',
                      user='USER', password='password')
    cur = con.cursor()

    exp_id = 832

    failed = False

    try:
        enz = load_model(cur, exp_id)
    except:
        failed = True
        print(traceback.format_exc())
    finally:
        cur.close()

    if not failed:
        print("Creating the Omex archive...")
        enz.create_archive()
        print("Archive created and finished.")
