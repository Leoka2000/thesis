from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from django.contrib import messages

def molecule_view(request, molecule_name):
    # Assuming you have a dictionary with molecule names and data (e.g., reagents_data or products_data)
    # Retrieve data from session or generate it
    reagents_data = request.session.get('reagents_data', {})
    products_data = request.session.get('products_data', {})

    # Find the molecule data from reagents_data or products_data
    molecule_data = reagents_data.get(molecule_name) or products_data.get(molecule_name)

    return render(request, 'pages/widgets/oxidoreductases_widgets/oxidoreductases_script.html', {
        'name': molecule_name,
        'smiles': molecule_data.get('smiles') if molecule_data else None,
        'image_base64': molecule_data.get('image_base64') if molecule_data else None,
    })
    
def oxidoreductases_view(request):
    # Define SMILES strings for reagents and products (Laccase reaction)
    laccase_reagents = {
        "4-Benzenediol": "C1=CC(=CC=C1O)O",
        "Oxygen": "O=O", 
    }
    laccase_products = {
        "Benzosemiquinone": "C1=CC(=CC=C1[O])O",
        "Water": "O",
    }

    # Define SMILES strings for reagents and products (Lignin Peroxidase reaction)
    lignin_peroxidase_reagents = {
        "Lignin Reactant": "COC1=C(C=C(C=C1)C(CO)(CO)OC2=CC=CC=C2OC)OC",
        "Hydrogen Peroxide": "OO",
    }
    lignin_peroxidase_products = {
        "3,4-Dimethoxybenzaldehyde": "COC1=CC=C(C=C1OC)C=O ",
        "2-Methoxyphenol": "COC1=CC=CC=C1O",
        "Glycolaldehyde": "O=CCO",
        "Water": "O",
    }

    # Define SMILES strings for reagents and products (Versatile Peroxidase reaction)
    versatile_peroxidase_reagents = {
        "Versatile Reactant": "COC1=C(C=C(C=C1)OC)CC(CO)OC1=CC=CC=C1OC",
        "Hydrogen Peroxide": "OO",
    }
    versatile_peroxidase_products = {
        "3,4-Dimethoxybenzaldehyde": "COC1=CC=C(C=C1OC)C=O",
        "2-Methoxyphenol": "COC1=CC=CC=C1O",
        "Glycolaldehyde": "O=CCO",
        "Water": "O",
    }

    # Define SMILES strings for reagents and products (Manganese Peroxidase reaction)
    manganese_peroxidase_reagents = {
        "Mn²⁺": "[Mn+2]",
        "Hydrogen Ions": "[H+].[H+]",
        "Hydrogen Peroxide": "OO",
    }
    manganese_peroxidase_products = {
        "Mn³⁺": "[Mn+3]",
        "Water": "O.O",
    }

    selected_enzyme = request.GET.get('specific_enzyme', '').lower()
    selected_chemical = request.GET.get('chemical_group', '').lower()

    def process_molecules(molecule_dict):
        processed_data = {}
        for name, smiles in molecule_dict.items():
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                # Add explicit hydrogen atoms to the molecule
                molecule_with_h = Chem.AddHs(molecule)

                # Convert the molecule to SMILES and generate the image
                smiles_string = Chem.MolToSmiles(molecule_with_h)

                # Generate the image with hydrogens explicitly shown
                img = Draw.MolToImage(molecule_with_h, kekulize=True, highlightAtoms=[], size=(500, 500))

                # Save image to a buffer in PNG format
                buffered = BytesIO()
                img.save(buffered, format="PNG")

                # Encode the image in base64
                img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                
                # Store the processed data
                processed_data[name] = {
                    "smiles": smiles_string,
                    "image_base64": img_base64,
                }
        return processed_data

    reaction_available = False
    reagents_data = {}
    products_data = {}

    # Check for reactions
    if request.method == "GET" and 'specific_enzyme' in request.GET and 'chemical_group' in request.GET:
        if selected_enzyme == "laccase" and selected_chemical == "carbonyl":
            reaction_available = True
            reagents_data = process_molecules(laccase_reagents)
            print(laccase_reagents)
            products_data = process_molecules(laccase_products)
        elif selected_enzyme == "lignin_peroxidase" and selected_chemical == "carbonyl":
            reaction_available = True
            reagents_data = process_molecules(lignin_peroxidase_reagents)
            products_data = process_molecules(lignin_peroxidase_products)
        elif selected_enzyme == "versatile_peroxidase" and selected_chemical == "carbonyl":
            reaction_available = True
            reagents_data = process_molecules(versatile_peroxidase_reagents)
            products_data = process_molecules(versatile_peroxidase_products)
        elif selected_enzyme == "manganese_peroxidase" and selected_chemical == "carbonyl":
            reaction_available = True
            reagents_data = process_molecules(manganese_peroxidase_reagents)
            products_data = process_molecules(manganese_peroxidase_products)

        # Save data to session
        request.session['reaction_available'] = reaction_available
        request.session['reagents_data'] = reagents_data
        request.session['products_data'] = products_data
        request.session['selected_enzyme'] = selected_enzyme

        # Messages
        if reaction_available:
            messages.success(request, "Reaction found successfully!")
        else:
            messages.error(request, "Reaction not available for selected enzyme or chemical group.")

    # If session has data, load it
    if 'reaction_available' in request.session and 'reagents_data' in request.session and 'products_data' in request.session:
        reaction_available = request.session['reaction_available']
        reagents_data = request.session['reagents_data']
        products_data = request.session['products_data']
        selected_enzyme = request.session.get('selected_enzyme', '')

    return render(request, "pages/oxidoreductases.html", {
        "reaction_available": reaction_available,
        "reagents_data": reagents_data,
        "products_data": products_data,
        "selected_enzyme": selected_enzyme,
    })
