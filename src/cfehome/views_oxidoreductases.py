from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from django.contrib import messages


def molecule_view(request, molecule_name):
     return render(request, 'pages/widgets/oxidoreductases_widgets/oxidoreductases_script.html', {'name': molecule_name})
def oxidoreductases_view(request):
    # Default SMILES for various reactions
    reactions = {
        'laccase': {
            'reagents': {
                "4-Benzenediol": "C1=CC(=CC=C1O)O",
                "Oxygen": "O=O",
            },
            'products': {
                "Benzosemiquinone": "C1=CC(=CC=C1[O])O",
                "Water": "O",
            }
        },
        'lignin_peroxidase': {
            'reagents': {
                "Lignin Reactant": "COC=1C=C(C=CC1OC)C(CC1=CC=C(C=C1)O)O",
                "Hydrogen Peroxide": "OO"
            },
            'products': {
                "3,4-Dimethoxybenzaldehyde": "COC1=CC=C(C=C1OC)C=O",
                "2-Methoxyphenol": "COC1=CC=CC=C1O",
                "Glycolaldehyde": "O=CCO",
                "Water": "O"
            }
        },
        'versatile_peroxidase': {
            'reagents': {
                "Versatile Reactant": "COC1=C(C=C(C=C1)OC)CC(CO)OC1=CC=CC=C1OC",
                "Hydrogen Peroxide": "OO"
            },
            'products': {
                "3,4-Dimethoxybenzaldehyde": "COC1=CC=C(C=C1OC)C=O",
                "2-Methoxyphenol": "COC1=CC=CC=C1O",
                "Glycolaldehyde": "O=CCO",
                "Water": "O"
            }
        },
        'manganese_peroxidase': {
            'reagents': {
                "Mn²⁺": "[Mn+2]",
                "Hydrogen Ions": "[H+].[H+]",
                "Hydrogen Peroxide": "OO"
            },
            'products': {
                "Mn³⁺": "[Mn+3]",
                "Water": "O.O"
            }
        }
    }

    selected_enzyme = request.GET.get('specific_enzyme', '').lower()
    selected_chemical = request.GET.get('chemical_group', '').lower()

    def generate_images(smiles_dict):
        images = {}
        for name, smi in smiles_dict.items():
            mol = Chem.MolFromSmiles(smi)
            img = Draw.MolToImage(mol)
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            encoded_image = base64.b64encode(buffer.getvalue()).decode("utf-8")
            images[name] = {"smiles": smi, "image_base64": encoded_image}
        return images

    reaction_available = False
    reagents_images = []
    products_images = []

    # Handle form submission and update session data
    if request.method == "GET" and 'specific_enzyme' in request.GET and 'chemical_group' in request.GET:
        if selected_enzyme in reactions and selected_chemical == "carbonyl":
            reaction_available = True
            reagents_images = generate_images(reactions[selected_enzyme]['reagents'])
            products_images = generate_images(reactions[selected_enzyme]['products'])

            # Store the results in session for later use
            request.session['reaction_available'] = reaction_available
            request.session['reagents_data'] = reagents_images
            request.session['products_data'] = products_images
            request.session['selected_enzyme'] = selected_enzyme

            # Display appropriate message based on reaction availability
            if reaction_available:
                messages.success(request, "Reaction found successfully!")
            else:
                messages.error(request, "Reaction not available for selected enzyme or chemical group.")
    
    # If session has data (from previous page load or refreshing), use that
    if 'reaction_available' in request.session and 'reagents_data' in request.session and 'products_data' in request.session:
        reaction_available = request.session['reaction_available']
        reagents_images = request.session['reagents_data']
        products_images = request.session['products_data']
        selected_enzyme = request.session.get('selected_enzyme', '')

    context = {
        "reaction_available": reaction_available,
        "reagents_data": reagents_images,
        "products_data": products_images,
        "selected_enzyme": selected_enzyme,  # for selecting enzyme in the form
    }

    return render(request, "pages/oxidoreductases.html", context)
