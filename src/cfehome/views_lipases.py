from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def lipases_view(request):
    # Define SMILES strings for reagents and products
    reagents = {
        "POP (1,3-Dipalmitoyl-2-oleylglycerol)": "CCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCC",
        "Stearic Acid": "CCCCCCCCCCCCCCC(=O)O",
        "Triolein": "O=C(OCC(OC(=O)CCCCCCC\\C=C/CCCCCCCC)COC(=O)CCCCCCC\\C=C/CCCCCCCC)CCCCCCC\\C=C/CCCCCCCC",
        "Palm Top Oil": "CCCCCCCCCCCCCCCC(=O)O"
    }
    products = {
        "P-OSt (1-palmitoyl-3-stearoyl-2-oleylglycerol)": "CCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCC",
        "St-OSt (1,3-distearoyl-2-oleylglycerol)": "CCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCC",
        "HMS": "C1C2(COP2O1)CSN",
        "Product": "C(C1CP(F)O1)F"
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
            # If the selected values match, process molecules and show success message
            specific_reagents = {
                "Palm Top Oil": reagents["Palm Top Oil"],
                "Triolein": reagents["Triolein"]
            }
            specific_products = {
                "HMS": products["HMS"],
                "Product": products["Product"]
            }
            reagents_data = process_molecules(specific_reagents)
            products_data = process_molecules(specific_products)
            messages.success(request, "Reaction successful! The reagents and products are displayed.")
        else:
            # If the selected values don't match, add an error message
            messages.error(request, "The selected options do not match the expected reaction.")

    return render(request, 'pages/lipases.html', {
        "reagents_data": reagents_data,
        "products_data": products_data
    })
