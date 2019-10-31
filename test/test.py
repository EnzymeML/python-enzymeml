import enzymeml.enzymeml as enzml
import libsbml as sbml

if __name__ == "__main__":
    enzml.DEBUG = False
    enzymeml = enzml.EnzymeML("test")

    enzymeml.add_creator("Mustermann", "Maxi", "uni", None)

    enzymeml.add(enzml.key.MAIN_META_CREATOR, {"family": "Mustermann", "given": "Max"})
    enzymeml.add(enzml.key.MAIN_META_DATES_CREATE, sbml.Date(2017, 6, 15))
    enzymeml.add(enzml.key.MAIN_META_DATES_MODIFY, None)
    unit = enzymeml.add(enzml.key.MAIN_UNIT,
                              {
                                  "name": "per_second", "units": [{"kind": sbml.UNIT_KIND_SECOND, "exponent": -1}]
                              })
    enzymeml.add(enzml.key.MAIN_UNIT_IS, ["resource:1", "resource:2"], unit)
    enzymeml.add(enzml.key.MAIN_UNIT_IS, "resource:3", unit)
    comp = enzymeml.add(enzml.key.MAIN_COMPARTMENT, {"name": "unidentified"})
    enzymeml.add(enzml.key.MAIN_COMPARTMENT_IS, "identifiers:go:0012", comp)
    s1 = enzymeml.add(enzml.key.MAIN_SPECIES,
                            {
                                "name": "n-pentanal", "compartment": comp, "type": enzml.ontology.SBO_SUBSTRATE,
                                "init_conc": 0.45, "units": unit
                            })
    enzymeml.add(enzml.key.MAIN_SPECIES_SPECIES, {"smiles": "CCCCC=O", "is": "chemi:005"}, s1)
    s2 = enzymeml.add(enzml.key.MAIN_SPECIES,
                            {
                                "name": "EctktA_I189G/H261V", "compartment": comp, "type": enzml.ontology.SBO_PRODUCT,
                                "init_conc": 0.45, "units": unit
                            })
    s3 = enzymeml.add(enzml.key.MAIN_SPECIES,
                            {
                                "name": "n-pentanal", "compartment": comp, "type": enzml.ontology.SBO_SUBSTRATE,
                                "init_conc": 0.22, "units": unit
                            })
    enzymeml.add(enzml.key.MAIN_SPECIES_PROTEIN,
                       {
                           "sequence": "Some fasta", "is": "pdb:006", "hasPart": ["pdb:0012", "pdb:0015"]
                       }, s2)
    r1 = enzymeml.add(enzml.key.MAIN_REACTION,
                            {
                                "name": "Some Reaction", "reactants": [{"id": s1, "stochiometry": 2}]
                            })

    enzymeml.add(enzml.key.MAIN_REACTION_CONDITION, {
        "temperature": (273.16, "kelvin"), "pressure": (101, "kPa"), "ph": 7.2
    }, r1)

    enzformat = enzml.EnzymeMLFormat()
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_TIME, "seconds"))
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_EMPTY))
    c1 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, s1, unit, "other", "The descr"))
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_EMPTY, 3, "Some comments"))
    c2 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, s1, unit))

    form1 = enzymeml.add(enzml.key.MAIN_DATA_FORMAT, enzformat)
    file1 = enzymeml.add(enzml.key.MAIN_DATA_FILE, {"format": form1, "file": "measures.csv"})
    measure1 = enzymeml.add(enzml.key.MAIN_DATA_MEASUREMENTS,
                                  {
                                      "file": file1, "start": 10, "stop": 55, "name": "Test measurement"
                                  })
    repl1 = enzml.EnzymeMLReplica(measure1, c1.replica)
    enzymeml.add(enzml.key.MAIN_REACTION_REPLICAS, repl1, r1)
    repl2 = enzml.EnzymeMLReplica(measure1, c2.replica)
    enzymeml.add(enzml.key.MAIN_REACTION_REPLICAS, repl2, r1)

    r2 = enzymeml.add(enzml.key.MAIN_REACTION,
                            {
                                "name": "Other Reaction", "reactants": [{"id": s1, "stochiometry": 2}]
                            })
    enzymeml.add(enzml.key.MAIN_REACTION_REPLICAS, enzml.EnzymeMLReplica(measure1, c1.replica), r2)
    enzymeml.add(enzml.key.MAIN_REACTION_REPLICAS, enzml.EnzymeMLReplica(measure1, c2.replica), r2)
    """

    reac1 = enzml.add_to_model(model, enzml.key.MAIN_REACTION, None, {"name": "reac"})

    enzformat = enzml.EnzymeMLFormat()
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_TIME, "seconds"))
    c1 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "S1", "mM", "R1"))
    c2 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "P1", "mM", "R1"))
    c3 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "S1", "mM", "R2"))
    c4 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "P1", "mM", "R2"))
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_EMPTY, 1))
    c5 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "S5", "mM", "R3"))
    c6 = enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_CONCENTRATION, "P5", "mM", "R3"))
    enzformat.add_column(enzml.create_column(enzml.COLUMN_TYPE_EMPTY))

    form1 = enzml.add_to_model(model, enzml.key.MAIN_DATA_FORMAT, None, enzformat)
    file1 = enzml.add_to_model(model, enzml.key.MAIN_DATA_FILE, None, {"format": form1, "file": "measures.csv"})
    measure1 = enzml.add_to_model(model, enzml.key.MAIN_DATA_MEASUREMENTS, None,
                                  {
                                      "file": file1, "start": 2, "stop": 4, "name": "Measure1"
                                  })
    measure2 = enzml.add_to_model(model, enzml.key.MAIN_DATA_MEASUREMENTS, None,
                                  {
                                      "file": file1, "start": 6, "stop": 7, "name": "Measure2"
                                  })# FIXME"""

    emod = enzymeml.create_model("testmodel")

    mr1 = emod.add(enzml.key.MODEL_REACTION,
                             {
                                 "name": "Some Reaction", "kineticlaw": "(S0*S1*Vm)/(Km*S2)",
                                 "parameters": [{"name": "Km", "value": 3.34, "units": unit, "stdev": 0.1},# TODO units
                                                {"name": "Vm", "value": 312.3,
                                                 "units": sbml.UNIT_KIND_SECOND,
                                                 "stdev": (1.2, sbml.UNIT_KIND_DIMENSIONLESS)}],
                                 "reactants": [{"id": s1}, {"id": s2}],
                                 "products": [{"id": s3, "stochiometry": 2}]
                             })

    emod.add(enzml.key.MODEL_REACTION_DATA, {r1: [repl1.id, repl2]}, mr1)
    enzymeml.create_archive()

    print("finished")
