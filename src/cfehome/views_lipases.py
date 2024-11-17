from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def lipases_view(request):
    # Define SMILES strings for reagents and products (First reaction)
    reagents = {
        "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",  # IUPAC: 2,3-bis(hexadecanoyloxy)propyl (Z)-octadec-9-enoate
        "Stearic Acid": "C(CCCCCCCCCCCCCCCCC)(=O)O",  # IUPAC: octadecanoic acid
    }
    products = {
        "P-OSt - Cocoa butter substute": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O",  # IUPAC: 2-hydroxy-1-(octadecanoyloxy)propyl (Z)-octadec-9-enoate
        "St-OSt": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
        "Palmitate": "C(CCCCCCCCCCCCCCC)(=O)[O-]",
    }

    # Define second set of SMILES strings for reagents and products (Second reaction)
    second_cbs_reaction_agents = {
        "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "St-St-St": "C(CCCCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",  # Additional reagent
    }
    second_cbs_reaction_products = {
        "P-OSt - Cocoa butter substute": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O",
        "St-O-St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-St": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-P": "C(CCCCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",  # Additional product
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

        # Check if selected values match the expected combination for the second reaction
        if (
            palm_oil == "palm_oil_midfraction"
            and organic_acid == "stearic_acid_3st"
            and catalyser == "1_3_specific_lipase"
        ):
            # Process molecules for the second reaction
            reagents_data = process_molecules(second_cbs_reaction_agents)
            products_data = process_molecules(second_cbs_reaction_products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")
        else:
            # Process molecules for the first reaction (default)
            reagents_data = process_molecules(reagents)
            products_data = process_molecules(products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")

    return render(request, 'pages/lipases.html', {
        "reagents_data": reagents_data,
        "products_data": products_data,
    })
