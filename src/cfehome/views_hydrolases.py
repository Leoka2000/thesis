from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

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

    # Check if the form was submitted by verifying the presence of query parameters
    form_submitted = 'added_value_molecule' in request.GET

    # Determine the molecule to display based on selected value in form
    value_added_molecule = request.GET.get('added_value_molecule', 'vanillin')
    added_molecule_image = None
    if form_submitted:
        if value_added_molecule == 'vanillin':
            added_molecule_image = generate_image(vanillin_smiles)
        elif value_added_molecule == '1_3_propadienol':
            added_molecule_image = generate_image(propadienol_smiles)
        elif value_added_molecule == 'cyclic_acetal':
            added_molecule_image = generate_image(cyclic_acetal_smiles)

        # Add a success or error message only if the form was submitted
        if added_molecule_image:
            messages.success(request, f'Successfully found the molecule: {value_added_molecule.replace("_", " ").title()}')
        else:
            messages.error(request, 'No matching molecule found for the selected value.')

    return render(request, "pages/hydrolases.html", {
        "pet_image": pet_image,
        "pet_polymer_image": pet_polymer_image,
        "added_molecule_image": added_molecule_image,
        "selected_molecule": value_added_molecule,
    })
