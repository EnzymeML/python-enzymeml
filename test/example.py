import libsbml as sbml
import enzymeml.enzymeml as enzymeml

doc = enzymeml.create_sbml_document()
model = doc.getModel()

enzymeml.add_to_model(model, enzymeml.key.MAIN_META_CREATOR, None, {"family": "Mustermann", "given": "Max"})
unit1 = enzymeml.add_to_model(model, enzymeml.key.MAIN_UNIT, None,
                              {
                                  "name": "mM", "units": [{"kind": sbml.UNIT_KIND_MOLE, "scale": -3},
                                                          {"kind": sbml.UNIT_KIND_LITRE, "exponent": -1}]
                              })
comp = enzymeml.add_to_model(model, enzymeml.key.MAIN_COMPARTMENT, None, {"name": "unidentified"})
s1 = enzymeml.add_to_model(model, enzymeml.key.MAIN_SPECIES, None,
                           {
                               "name": "n-pentanal", "compartment": comp, "type": enzymeml.ontology.SBO_SUBSTRATE,
                               "init_conc": 0.45, "units": unit1
                           })
enzymeml.add_to_model(model, enzymeml.key.MAIN_SPECIES_SPECIES, s1, {"smiles": "CCCCC=O",
                                                                     "is": "http://identifiers.org/chebi/CHEBI:84069"})
r1 = enzymeml.add_to_model(model, enzymeml.key.MAIN_REACTION, None,
                           {
                               "name": "Some Reaction", "reactants": [enzymeml.EnzymeMLReactionComponent(s1, 2)]
                           })
enzymeml.add_to_model(model, enzymeml.key.MAIN_REACTION_CONDITION, r1, {"temperature": (125, sbml.UNIT_KIND_KELVIN),
                                                                        "ph": 7.2})

sbml.writeSBMLToFile(doc, "ex.xml")
