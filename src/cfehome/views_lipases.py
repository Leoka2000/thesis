# views.py
from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def lipases_view(request):
    # Define SMILES strings for reagents and products
    reagents = {
        "POP (1,3-Dipalmitoyl-2-oleylglycerol)": "C(CCCCCCCCCCCCCCC)(=O)OCC(OCCCCCCCC\C=C/CCCCCCCC)COC(CCCCCCCCCCCCCCC)=O",
        "Stearic Acid": "C(CCCCCCCCCCCCCCCCC)(=O)O",
        "palm_top_fraction (tripalmitin)": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "Triolein": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O",
    }
    products = {
        "P-OSt - CBS - (2-stearoyl-1,3-dihydroxypropane)": "C(CCCCCCCCCCCCCCCCC)(=O)OC(CO)CO",
        "St-OSt (1,3-distearoyl-2-oleylglycerol)": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-St": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
        "P-St-P": "C(CCCCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
        "HMS (P-O-O)": "C(CCCCCCCCCCCCCCC)(=O)OC[C@H](COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O",
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
                    "image_base64": img_base64
                }
        return processed_data

    # Default empty data for rendering
    reagents_data = {}
    products_data = {}

    # Process request parameters if they exist
    if request.method == "GET":
        palm_oil = request.GET.get("palm_oil_mildfraction")
        organic_acid = request.GET.get("organic_fatty_acid")
        catalyser = request.GET.get("catalyser")

        # Check if selected values match the expected combination
        if (
            palm_oil == "palm_top_oil"
            and organic_acid == "triolein"
            and catalyser == "1_3_specific_lipase"
        ):
            # Process specific molecules for this reaction
            specific_reagents = {
                "palm_top_fraction (tripalmitin)": reagents["palm_top_fraction (tripalmitin)"],
                "Triolein": reagents["Triolein"]
            }
            specific_products = {
                "HMS (P-O-O)": products["HMS (P-O-O)"],
                "P-O-P": products["P-O-P"]
            }
            reagents_data = process_molecules(specific_reagents)
            products_data = process_molecules(specific_products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")
        elif (
            palm_oil == "palm_oil_mildfraction"
            and organic_acid == "stearic_acid"
            and catalyser == "1_3_specific_lipase"
        ):
            # Process default reaction
            reagents_data = process_molecules(reagents)
            products_data = process_molecules(products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")
        else:
            # If the selected values don't match, add an error message
            messages.error(request, "The selected options do not match the expected reaction.")

    return render(request, 'pages/lipases.html', {
        "reagents_data": reagents_data,
        "products_data": products_data
    })

