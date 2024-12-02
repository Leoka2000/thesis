
from django.contrib import admin
from django.urls import path

from . import views
from . import views_oxidoreductases
from . import views_lipases
from . import views_hydrolases

urlpatterns = [
    path('', views.home_view, name='home'),
    path('lipases', views_lipases.lipases_view, name='lipases'),
    path('hydrolases', views_hydrolases.hydrolases_view, name='hydrolases'),
    path('transferases', views.transferases_view, name='transferases'),
    path('oxidoreductases', views_oxidoreductases.oxidoreductases_view, name='oxidoreductases'),
    path('molecule/<str:molecule_name>/', views_lipases.molecule_detail_view, name='molecule_detail'),
  
    path('admin/', admin.site.urls),
    
]

