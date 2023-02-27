'''
Gestion des éléments tiges
'''

def extract_stems_from_MTG(MTG, entity_id):
    stems=[]
    geom = MTG.property('geometry')
    for vid in geom.keys():
        if MTG.class_name(vid) == 'StemElement':
            stems.append((vid, entity_id))
    
    return stems

def manage_stems_for_ratp(stems_id, matching_ids, ratp_parameters) :
    if len(stems_id) > len(matching_ids) :
        raise ValueError("Too many stems elements ")
    
    current_number_of_entities = max([v[1] for v in matching_ids.values()]) + 1
    mem_mu=[]
    for couple in stems_id:
        # on recherche la shape correspondante dans le dict des id
        glob_id = [k for k,v in matching_ids.items() if tuple(v) == couple][0]

        # change le numéro de l'entité
        old_values = matching_ids[glob_id]
        matching_ids[glob_id] = [old_values[0],
                                    old_values[1] + current_number_of_entities]

        # ajoute des propriétés optiques et mu sur la nouvelle entité tige
        if couple[1] not in mem_mu:
            ratp_parameters["mu"].append(ratp_parameters["mu"][couple[1]])
            element = ratp_parameters["reflectance coefficients"][couple[1]]
            ratp_parameters["reflectance coefficients"].append(element)
            mem_mu.append(couple[1])