from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

@csrf_exempt
def submit_molecule(request):
    if request.method == 'POST':
        molecule_name = request.POST.get('molecule_name')
        print(f"Received molecule: {molecule_name}")
        # Store the molecule name in the session
        request.session['molecule_name'] = molecule_name
        return JsonResponse({'status': 'success', 'molecule': molecule_name})
    return JsonResponse({'status': 'error', 'message': 'Invalid request'}, status=400)

def lipases_view(request):
    molecule_name = request.session.get('molecule_name')
    # Define SMILES strings for reagents and products (First reaction)
    first_cbs_reaction_reagents = {
        "palm_oil_mid_fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "stearic_acid": "C(CCCCCCCCCCCCCCCCC)(=O)O",
    }
    first_cbs_reaction_products = {
        "P_OSt": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O", #Cocoa butter substitute
        "St_O_St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
        "palmitate": "C(CCCCCCCCCCCCCCC)(=O)[O-]",
    }

    # Second reaction
    second_cbs_reaction_reagents = {
        "palm_oil_mid_fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
        "St_St_St": "C(CCCCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
    }
    second_cbs_reaction_products = {
        "P_OSt": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(C(C)O)OC(CCCCCCCCCCCCCCCCC)=O",  #Cocoa butter substitute
        "St_O_St": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCCCCO))COC(CCCCCCCCCCCCCCCCC)=O",
        "P_St_St": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O",
        "P_St_P": "C(CCCCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
    }

    # HMS Triolein Reaction
    first_hms_triolein_reagents = {
        "palm_oil_top_fraction": "C(CCCCCCCCCCCCCCC)(=O)O",
        "triolein": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O",
    }
    first_hms_triolein_products = {
        "O_P_O": "C(CCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCC\C=C/CCCCCCCC)=O)COC(CCCCCCC\C=C/CCCCCCCC)=O", #HMS
        "P_O_P": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OC[C@@H](COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O",
    }

    second_hms_oleic_acid_reagents = {
        "palm_oil_top_fraction": "C(CCCCCCCCCCCCCCC)(=O)O",
        "oleic_acid_3O": "C(CCCCCCCC=CCCCCCCCC)(=O)OCC(COC(CCCCCCCC=CCCCCCCCC)=O)OC(CCCCCCCC=CCCCCCCCC)=O",
    }
    second_hms_oleic_acid_products = {
        "O_P_O": "C(CCCCCCCCCCCCCCC)(=O)OC(COC(CCCCCCC\\C=C/CCCCCCCC)=O)COC(CCCCCCC\\C=C/CCCCCCCC)=O", #HMS
        "P_P_O": "C(CCCCCCCCCCCCCCC)(=O)OCC(OC(CCCCCCC\\C=C/CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
        "palmitate": "C(CCCCCCCCCCCCCCC)(=O)[O-]",
    }


    def process_molecules(molecule_dict):
        processed_data = {}
        for name, smiles in molecule_dict.items():
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                smiles_string = Chem.MolToSmiles(molecule)
                img = Draw.MolToImage(molecule, size=(500, 200))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
                processed_data[name] = {
                    "smiles": smiles_string,
                    "image_base64": img_base64,
                }
        return processed_data

    # creating the "state", that is conteyt thatb we will pass
    reagents_data = {}
    products_data = {}
    show_reaction_image_top_oil_triolein = False
    show_reaction_image_midfraction_stearic_acid = False
    form_submitted = False

    
    # Pprocess mages if get request s succcessfull
    if request.method == "GET" and "palm_oil_midfraction" in request.GET:
        form_submitted = True
        palm_oil = request.GET.get("palm_oil_midfraction")
        organic_acid = request.GET.get("organic_fatty_acid")
        catalyser = request.GET.get("catalyser")

        # Check for matching reaction based on user input
        if (
            palm_oil == "palm_top_oil"
            and organic_acid == "triolein"
            and catalyser == "1_3_specific_lipase"
        ):
            show_reaction_image_top_oil_triolein = True
            reagents_data = process_molecules(first_hms_triolein_reagents)
            products_data = process_molecules(first_hms_triolein_products)
            messages.success(request, "HMS-Triolein reaction successful! The reagents and products are displayed.")
        elif (
            palm_oil == "palm_oil_midfraction"
            and organic_acid == "stearic_acid_3st"
            and catalyser == "1_3_specific_lipase"
        ):
            show_reaction_image_midfraction_stearic_acid = True
            reagents_data = process_molecules(second_cbs_reaction_reagents)
            products_data = process_molecules(second_cbs_reaction_products)
            messages.success(request, "Second P-OStReaction successful! The reagents and products are displayed.")
        elif (
            palm_oil == "palm_top_oil"
            and organic_acid == "oleic_acid"
            and catalyser == "1_3_specific_lipase"
        ):
            reagents_data = process_molecules(second_hms_oleic_acid_reagents)
            products_data = process_molecules(second_hms_oleic_acid_products)
            messages.success(request, "HMS-Oleic acid reaction successful! The reagents and products are displayed.")
        elif (
            palm_oil == "palm_oil_midfraction"
            and organic_acid == "stearic_acid"
            and catalyser == "1_3_specific_lipase"
        ):
            reagents_data = process_molecules(first_cbs_reaction_reagents)
            products_data = process_molecules(first_cbs_reaction_products)
            messages.success(request, "P-OSt reaction successful! The reagents and products are displayed.")
        else:
            messages.error(request, "Invalid input. Please check your selections.")

            
    return render(request, 'pages/lipases.html', {
        "reagents_data": reagents_data,
        "products_data": products_data,
        "show_reaction_image_top_oil_triolein": show_reaction_image_top_oil_triolein,
        "show_reaction_image_midfraction_stearic_acid": show_reaction_image_midfraction_stearic_acid,
         "form_submitted": form_submitted,
          "molecule_name": molecule_name,
   
       
        
    })




  # Define SMILES strings for reagents and products
    # reagents = {
    #     #POP= palm oil mid fraction
    #     "POP Palm oil mid fraction": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O", #iupac 2,3-bis(hexadecanoyloxy)propyl (Z)-octadec-9-enoate
    #     "Stearic Acid": "C(CCCCCCCCCCCCCCCCC)(=O)O", #IUPAC octadecanoic acid
    #     "palm_top_fraction": "C(CCCCCCCCCCCCCCC)(=O)O", #IUPAC Hexadecanoic acid
    #     #acima fiz qnd tava doente e arrumei
    #     "Triolein": "C(CCCCCCC\C=C/CCCCCCCC)(=O)OCC(COC(CCCCCCC\C=C/CCCCCCCC)=O)OC(CCCCCCC\C=C/CCCCCCCC)=O", #IUPAC 2,3-bis((Z)-octadec-9-enoyloxy)propyl (Z)-octadec-9-enoate
    #     "tripalmitin_tryglyceryde_3_oxygens": "C(CCCCCCCCCCCCCCC)(=O)OCC(COC(CCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCC)=O", #iupac 2,3-bis(hexadecanoyloxy)propyl hexadecanoate
    #     "simples_olete_acid": "C(CCCCCCC\C=C/CCCCCCCC)(=O)O" #(9Z)-octadec-9-enoic acid
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


    #iupac name for oleic acid with three oleates: Propane-1,2,3-triyl tris(9-octadecenoate)
    #iupac for P-P-O: 1,3-Dipalmitoyl-2-oleoylglycerol