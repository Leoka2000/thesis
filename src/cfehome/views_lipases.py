from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def lipases_view(request):
    # Define SMILES strings for reagents and products (First reaction)
    reagents = {
        "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "Stearic Acid": "C(CCCCCCCCCCCCCCCCC)(=O)O",
    }
    products = {
        "P-OSt - Cocoa butter substitute": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O",
        "St-O-St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
        "Palmitate": "C(CCCCCCCCCCCCCCC)(=O)[O-]",
    }

    # Second reaction
    second_cbs_reaction_agents = {
        "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "St-St-St": "C(CCCCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
    }
    second_cbs_reaction_products = {
        "P-OSt - Cocoa butter substitute": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O",
        "St-O-St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCCO)COC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-St": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-P": "C(CCCCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
    }

    # HMS Triolein Reaction
    hms_triolein_reagents = {
        "palm_top_fraction (tripalmitin)": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "Triolein": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O",
    }
    hms_triolein_products = {
        "HMS (O-P-O)": "C(CCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCC\C=C/CCCCCCCC)=O)COC(CCCCCCC\C=C/CCCCCCCC)=O",
        "P-O-P": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC[C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
    }

    def process_molecules(molecule_dict):
        processed_data = {}
        for name, smiles in molecule_dict.items():
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                smiles_string = Chem.MolToSmiles(molecule)
                img = Draw.MolToImage(molecule, size=(300, 300))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                processed_data[name] = {
                    "smiles": smiles_string,
                    "image_base64": img_base64,
                }
        return processed_data

    # Default empty data for rendering
    reagents_data = {}
    products_data = {}

    # Process request parameters if they exist
    if request.method == "GET":
        palm_oil = request.GET.get("palm_oil_midfraction")
        organic_acid = request.GET.get("organic_fatty_acid")
        catalyser = request.GET.get("catalyser")

        if (
            palm_oil == "palm_top_oil"
            and organic_acid == "triolein"
            and catalyser == "1_3_specific_lipase"
        ):
            reagents_data = process_molecules(hms_triolein_reagents)
            products_data = process_molecules(hms_triolein_products)
            messages.success(request, "HMS-Triolein reaction successful! The reagents and products are displayed.")
        elif (
            palm_oil == "palm_oil_midfraction"
            and organic_acid == "stearic_acid_3st"
            and catalyser == "1_3_specific_lipase"
        ):
            reagents_data = process_molecules(second_cbs_reaction_agents)
            products_data = process_molecules(second_cbs_reaction_products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")
        else:
            reagents_data = process_molecules(reagents)
            products_data = process_molecules(products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")

    return render(request, 'pages/lipases.html', {
        "reagents_data": reagents_data,
        "products_data": products_data,
    })



  # Define SMILES strings for reagents and products
    # reagents = {
    #     #POP= palm oil mid fraction
    #     "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O", #iupac 2,3-bis(hexadecanoyloxy)propyl (Z)-octadec-9-enoate
    #     "Stearic Acid": "C(CCCCCCCCCCCCCCCCC)(=O)O", #IUPAC octadecanoic acid
    #     "palm_top_fraction (tripalmitin)": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O", #IUPAC 2,3-bis(hexadecanoyloxy)propyl hexadecanoate
    #     #acima fiz qnd tava doente e arrumei
    #     "Triolein": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O", #IUPAC 2,3-bis((Z)-octadec-9-enoyloxy)propyl (Z)-octadec-9-enoate
    #     "oleic-acid": "C(CCCCCCC\C=C/CCCCCCCC)(=O)O", #iupac (Z)-octadec-9-enoic acid
    # }
    # products = {
    #     "P-OSt - CBS - (2-stearoyl-1,3-dihydroxypropane)": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O", #IUPAC: 2-hydroxy-1-(octadecanoyloxy)propyl (Z)-octadec-9-enoate
    #     "St-OSt (1,3-distearoyl-2-oleylglycerol)": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(COC(CCCCCCCCCCCCCCCCC)=O)O)OC(CCCCCCCCCCCCCCCCC)=O", # IUPAC 2-hydroxy-1,3-bis(octadecanoyloxy)propyl (Z)-octadec-9-enoate
    #     "St-O-St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",#2-((Z)-octadec-9-enoyloxy)-1,3-bis(octadecanoyloxy)propane
    #      "P": "C(CCCCCCCCCCCCCCC)(=O)[O-]", #hexadecanoate 
    #     "P-St-St": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O", #IUPAC: 2,3-bis(octadecanoyloxy)propyl hexadecanoate
    #     "P-St-P": "C(CCCCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O", #IUPAC: 2-(octadecanoyloxy)-1,3-bis(hexadecanoyloxy)propane
    #     "HMS (O-P-O)": "C(CCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCC\C=C/CCCCCCCC)=O)COC(CCCCCCC\C=C/CCCCCCCC)=O", #iupac 2-(hexadecanoyloxy)-1,3-bis((Z)-octadec-9-enoyloxy)propane
    #     "P-O-P": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC[C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O", #IUPAC: (2R)-3-(hexadecanoyloxy)-2-(hexadecanoyloxy)propyl (9Z)-octadec-9-enoate
    #     "St-St-St": "C(CCCCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O"  #a trygliceride with three stereate group IUPAC: 2,3-bis(octadecanoyloxy)propyl octadecanoate
    # }