from cathub.cathubsql import CathubSQL
db = CathubSQL()

pub_id = "TangNitrate2024"
dataframe = db.get_dataframe(pub_id=pub_id, include_atoms=False)
print(dataframe)
dataframe.to_pickle(pub_id + '.pickle')

import pandas
print(dataframe.columns)
print(dataframe[['chemical_composition', 'surface_composition','facet', 'equation', 'reaction_energy']].to_markdown())
