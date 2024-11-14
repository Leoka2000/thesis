from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from django.contrib import messages

def oxidoreductases_view(request):
    # Default SMILES for an existing reaction
    reagents_smiles = ["C1=CC(=CC=C1O)O", "O=O"]  # 4-benzenediol (1,2-benzenediol) and O₂
    products_smiles = ["C1=CC(=CC=C1[O])O", "O"]   # Benzosemiquinone and H₂O
    
    # Additional SMILES for the lignin reaction
    lignin_reaction_reagents = [
        "COC1=C(C=C(C=C1)OC)CC(CO)O",  # 1-(3,4-dimethoxyphenyl)-2-(2-methoxyphenoxy)propane-1,3-diol
        "OO"  # Hydrogen peroxide (H₂O₂)
    ]
    lignin_reaction_products = [
        "COC1=CC=C(C=C1OC)C=O",        # 3,4-dimethoxybenzaldehyde
        "COC1=CC=CC=C1O",              # 2-methoxyphenol
        "O=CCO",                       # Glycolaldehyde
        "O"                            # Water (H₂O)
    ]
    
    # SMILES for versatile peroxidase reaction
    versatile_reaction_reagents = [
        "COC1=C(C=C(C=C1)OC)CC(CO)OC1=CC=CC=C1OC",  # 1-(3,4-dimethoxyphenyl)-2-(2-methoxyphenoxy)propane-1,3-diol
        "OO"  # Hydrogen peroxide (H₂O₂)
    ]
    versatile_reaction_products = [
        "COC1=CC=C(C=C1OC)C=O",        # 3,4-dimethoxybenzaldehyde
        "COC1=CC=CC=C1O",              # 2-methoxyphenol
        "O=CCO",                       # Glycolaldehyde
        "O"                            # Water (H₂O)
    ]
    
    # New SMILES for manganese peroxidase reaction
    manganese_reaction_reagents = [
        "[Mn+2]",      # Mn²⁺
        "[H+].[H+]",   # 2H⁺
        "OO"           # H₂O₂ (hydrogen peroxide)
    ]
    manganese_reaction_products = [
        "[Mn+3]",      # Mn³⁺
        "O.O"          # 2H₂O (water molecules)
    ]

    # Set default for no reaction available
    reaction_available = False
    reagents_images = []
    products_images = []

    # Check the selected enzyme and chemical group
    selected_enzyme = request.GET.get('specific_enzyme', '').lower()
    selected_chemical = request.GET.get('chemical_group', '').lower()
    
    # Function to generate images from SMILES
    def generate_images(smiles_list):
        images = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            img = Draw.MolToImage(mol)
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            encoded_image = base64.b64encode(buffer.getvalue()).decode("utf-8")
            images.append(encoded_image)
        return images

    # Check conditions for different reactions
    if selected_enzyme == "laccase" and selected_chemical == "carbonyl":
        reaction_available = True
        reagents_images = generate_images(reagents_smiles)
        products_images = generate_images(products_smiles)
    
    elif selected_enzyme == "lignin_peroxidase" and selected_chemical == "carbonyl":
        reaction_available = True
        reagents_images = generate_images(lignin_reaction_reagents)
        products_images = generate_images(lignin_reaction_products)
        
    elif selected_enzyme == "versatile_peroxidase" and selected_chemical == "carbonyl":
        reaction_available = True
        reagents_images = generate_images(versatile_reaction_reagents)
        products_images = generate_images(versatile_reaction_products)
        
    elif selected_enzyme == "maganese_peroxidade" and selected_chemical == "carbonyl":  # Added new condition for manganese peroxidase
        reaction_available = True
        reagents_images = generate_images(manganese_reaction_reagents)
        products_images = generate_images(manganese_reaction_products)
    if reaction_available:
        messages.success(request, "Reaction found successfully!")
    else:
        messages.error(request, "Reaction not available for selected enzyme or chemical group.")
    context = {
        "reaction_available": reaction_available,
        "reagents_images": reagents_images,
        "products_images": products_images,
    }
    
    return render(request, "pages/oxidoreductases.html", context)