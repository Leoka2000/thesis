from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def hydrolases_view(request):
    # Define the SMILES strings for each reaction
    reactions = {
        "first_reaction": {
            "pet_monomer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCC',
            "pet_polymer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCCC(=O)C2=CC=CC=C2C(=O)OCCC(=O)C3=CC=CC=C3C(=O)OCC',
            "vanillin": 'COC1=CC=C(C=C1)C=O',
        },
        "second_reaction": {
            "pet_monomer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCC',
            "pet_polymer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCCC(=O)C2=CC=CC=C2C(=O)OCCC(=O)C3=CC=CC=C3C(=O)OCC',
            "cyclic_acetal": 'O1CCOC1',
        },
        "third_reaction": {
            "pet_monomer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCC',
            "pet_polymer": 'O=C(OCC)C1=CC=CC=C1C(=O)OCCC(=O)C2=CC=CC=C2C(=O)OCCC(=O)C3=CC=CC=C3C(=O)OCC',
            "propadienol": 'C=CCO',
        }
    }

    # Generate the images for molecules
    def generate_image(smiles, size=(200, 200)):
        mol = Chem.MolFromSmiles(smiles)
        buffer = BytesIO()
        image = Draw.MolToImage(mol, size=size)
        image.save(buffer, format="PNG")
        return base64.b64encode(buffer.getvalue()).decode("utf-8")

    # Store generated images in a structure
    reaction_images = {}
    for reaction_name, molecules in reactions.items():
        reaction_images[reaction_name] = {
            name: generate_image(smiles) for name, smiles in molecules.items()
        }

    # Check if the form was submitted by verifying the presence of query parameters
    form_submitted = 'added_value_molecule' in request.GET

    # Determine the molecule to display based on selected value in form
    value_added_molecule = request.GET.get('added_value_molecule', 'vanillin')
    added_molecule_image = None

    if form_submitted:
        # Iterate over all reactions to find the selected molecule
        for molecules in reactions.values():
            if value_added_molecule in molecules:
                added_molecule_image = generate_image(molecules[value_added_molecule])
                break

        # Add a success or error message
        if added_molecule_image:
            messages.success(request, f'Successfully found the molecule: {value_added_molecule.replace("_", " ").title()}')
        else:
            messages.error(request, 'No matching molecule found for the selected value.')

    return render(request, "pages/hydrolases.html", {
        "reaction_images": reaction_images,
        "added_molecule_image": added_molecule_image,
        "selected_molecule": value_added_molecule,
    })
