from django.shortcuts import render
from django.contrib import messages
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def home_view(request):
    
    return render(request, "pages/home.html", {
        
    })



def transferases_view(request):
    return render(request, "pages/transferases.html", {})
