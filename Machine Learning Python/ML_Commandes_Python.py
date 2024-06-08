# Affectation de variables
x , y = 1 , 2.3

# Portes Logiques
print(False & True) # AND
print(False | True) # OR
print(False ^ True) # XOR

# Fonction courte
fonction_courte = lambda x , y : x**2 + y
res = fonction_courte(x = 3 , y = 2)
print(res)

## SEQUENCES
# Tuple - Structure de données :
# - à taille fixe
# - à données constantes
# => MEMOIRE ECONOMISEE SI TRAITEMENT BCP DE DATA
tuple_1 = (1 , 2 , 3 , 4)

# Slicing
Tableau = [1 , 2 , 3 , 4 , 5]
ind_deb = 2
ind_fin = 3
pas = 1
print(Tableau[ind_deb : ind_fin : pas])

# Insertion
ind_insert = 2
Tableau.insert( ind_insert , 10)

# Tri inversé
Tableau.sort(reverse = True)

# Index et Enumerate
for index , valeur in enumerate(Tableau) :
    print(index , valeur)
    
