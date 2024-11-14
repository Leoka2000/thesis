from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def home_view(request):
    # Define the SMILES strings for glucose monomers
    glucose_monomer_smiles = 'C(C1C(C(C(O1)O)O)O)O'
    glucose_monomer = Chem.MolFromSmiles(glucose_monomer_smiles)

    # Generate the image for each monomer and the combined form
    images = []
    for mol in [glucose_monomer, glucose_monomer]:  # To illustrate the two starting glucose monomers
        buffer = BytesIO()
        img = Draw.MolToImage(mol, size=(200, 200))
        img.save(buffer, format="PNG")
        images.append(base64.b64encode(buffer.getvalue()).decode("utf-8"))

    # Illustrate the polymerization (simplified view of two glucose monomers linked)
    # Create a placeholder molecule representing the linked form for visualization purposes
    polymer_smiles = 'C(C1C(C(C(O1)O)O)O)OC(C2C(C(C(O2)O)O)O)O'
    polymer_molecule = Chem.MolFromSmiles(polymer_smiles)
    polymer_image = Draw.MolToImage(polymer_molecule, size=(300, 300))

    # Convert to base64 for the template
    buffer = BytesIO()
    polymer_image.save(buffer, format="PNG")
    polymer_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return render(request, "pages/home.html", {
        "monomer_images": images,
        "polymer_image": polymer_img_str
    })



def hydrolases_view(request):
    # Define the SMILES strings
    pet_monomer_smiles = 'O=C(OCC)C1=CC=CC=C1C(=O)OCC'
    pet_polymer_smiles = 'O=C(OCC)C1=CC=CC=C1C(=O)OCCC(=O)C2=CC=CC=C2C(=O)OCCC(=O)C3=CC=CC=C3C(=O)OCC'
    vanillin_smiles = 'COC1=CC=C(C=C1)C=O'
    cyclic_acetal_smiles = 'O1CCOC1'
    propadienol_smiles = 'C=CCO'  # SMILES for 1,3-propadienol

    # Generate the images for PET polymer and PET monomer
    def generate_image(smiles, size=(200, 200)):
        mol = Chem.MolFromSmiles(smiles)
        buffer = BytesIO()
        image = Draw.MolToImage(mol, size=size)
        image.save(buffer, format="PNG")
        return base64.b64encode(buffer.getvalue()).decode("utf-8")

    pet_image = generate_image(pet_monomer_smiles, size=(200, 200))
    pet_polymer_image = generate_image(pet_polymer_smiles, size=(400, 400))

    # Determine the molecule to display based on selected value in form
    value_added_molecule = request.GET.get('added_value_molecule', 'vanillin')
    if value_added_molecule == 'vanillin':
        added_molecule_image = generate_image(vanillin_smiles)
    elif value_added_molecule == '1_3_propadienol':
        added_molecule_image = generate_image(propadienol_smiles)
    elif value_added_molecule == 'cyclic_acetal':
        added_molecule_image = generate_image(cyclic_acetal_smiles)
    else:
        added_molecule_image = None

    # Add a success or error message
    if added_molecule_image:
        messages.success(request, f'Successfully found the molecule: {value_added_molecule.replace("_", " ").title()}')
    else:
        messages.error(request, 'No matching molecule found for the selected value.')

    return render(request, "pages/hydrolases.html", {
        "pet_image": pet_image,
        "pet_polymer_image": pet_polymer_image,
        "added_molecule_image": added_molecule_image,
        "selected_molecule": value_added_molecule
    })


def transferases_view(request):
    return render(request, "pages/transferases.html", {})
