from typing import List, Optional

# complexité algorithme Simplexe

# Complexité Temporelle

#     Initialisation :
#         La création de la matrice des coefficients du bloc des contraintes et l'ajout des colonnes pour les variables d'écart ont une complexité de O(m×n)O(m×n), où mm est le nombre de contraintes et nn est le nombre de variables.
#         Le remplissage des colonnes des variables d'écart a une complexité de O(m)O(m).

#     Calcul des vecteurs ZjZj et Cj−ZjCj−Zj :
#         Le calcul de ZjZj et Cj−ZjCj−Zj a une complexité de O(n)O(n).

#     Sélection de la variable entrante :
#         La recherche de la valeur maximale dans Cj−ZjCj−Zj a une complexité de O(n)O(n).

#     Calcul des ratios et sélection de la variable sortante :
#         Le calcul des ratios et la recherche du ratio minimal ont une complexité de O(m)O(m).

#     Mise à jour de la matrice et des vecteurs :
#         La mise à jour de la ligne pivot et des autres lignes de la matrice a une complexité de O(m×n)O(m×n).
#         La mise à jour du vecteur QQ a une complexité de O(m)O(m).

#     Itérations :
#         Le nombre d'itérations nécessaires pour converger vers la solution optimale peut varier. Dans le pire des cas, la méthode du Simplexe peut nécessiter un nombre exponentiel d'itérations (complexité exponentielle), mais en pratique, elle converge généralement en un nombre polynomial d'itérations.

# Complexité Spatiale
#     La complexité spatiale est principalement déterminée par la taille de la matrice des contraintes et des vecteurs utilisés. La matrice des contraintes a une taille de O(m×(n+m))O(m×(n+m)) (en tenant compte des variables d'écart), et les vecteurs ZjZj, Cj−ZjCj−Zj, coef_pivotcoef_pivot, vdbvdb, et ratioratio ont une taille de O(n+m)O(n+m).

# Conclusion
#     Complexité temporelle : Dans le pire des cas, la complexité temporelle de la méthode du Simplexe est exponentielle, mais elle est souvent polynomiale.
#     Complexité spatiale : La complexité spatiale est de O(m×(n+m))O(m×(n+m)), ce qui est gérable pour des problèmes de taille modérée.
# =================================================================================================

# Fonction pour résoudre un problème d'optimisation linéaire Simplexe
# coef: Coefficients de la fonction objectif
# matrice_sous_contrainte: Matrice des contraintes
# vecteur: Côtés droits des inégalités
# type_inegalites: Type des inégalités
# maximisation: Maximiser ou minimiser (Sachant qu'on ne vas pas implémenter la minimisation)
def simplexe(
    coef: List[float],
    matrice_sous_contrainte: List[List[float]],
    vecteurs_Q,
    type_inegalite,
    maximisation=True,
):
    Z = 0

    # Création de la matrice des coefficients du bloc des contraintes
    matrice_coef_bloc_contrainte = matrice_sous_contrainte.copy()

    # Ajout des colonnes pour les variables d'écart
    for row in matrice_coef_bloc_contrainte:
        for _ in range(len(matrice_sous_contrainte)):
            row.append(0)

    # Remplissage des colonnes des variables d'écart
    for row_index in range(len(matrice_sous_contrainte)):
        if type_inegalite[row_index] == "<=":
            # On ajoute 1 pour les variables d'écart si c'est une inégalité <=
            matrice_coef_bloc_contrainte[row_index][
                row_index + len(matrice_sous_contrainte)
            ] = 1
        else:
            # On ajoute -1 pour les variables d'écart si c'est une inégalité >=
            matrice_coef_bloc_contrainte[row_index][
                row_index + len(matrice_sous_contrainte)
            ] = -1

    # Ajout des colonnes pour les variables d'écart
    for _ in range(len(matrice_sous_contrainte)):
        coef.append(0)

    # Création du vecteur Zj
    Zj = [0 for _ in range(len(coef))]

    # Création du vecteur Cj - Zj
    difference_Cj_Zj = [0 for _ in range(len(coef))]

    # Création du vecteur des coefficients pivot
    coef_pivot = [0 for _ in range(len(matrice_sous_contrainte))]

    for i in range(len(coef)):
        difference_Cj_Zj[i] = coef[i] - Zj[i]

    valeur_entrante = max(difference_Cj_Zj)
    colonne_valeur_entrante = difference_Cj_Zj.index(valeur_entrante)

    # Initialisation des variables dans base
    vdb = ["" for _ in range(len(coef_pivot))]
    for i in range(len(vdb)):
        vdb[i] = "x" + str(i + len(coef_pivot) + 1)

    # Initialisation des ratios
    ratio = [0 for _ in range(len(matrice_coef_bloc_contrainte))]

    # ============== Affichage du tableau Simplexe Initialisation =================
    display_tab(
        coef,
        matrice_coef_bloc_contrainte,
        vecteurs_Q,
        coef_pivot,
        Z,
        Zj,
        difference_Cj_Zj,
        vdb,
        ratio,
    )
    # ===================== Fin de l'affichage =====================

    # Calcul de la valeur sortante
    min_value = vecteurs_Q[0] / matrice_coef_bloc_contrainte[0][colonne_valeur_entrante]
    index_min_value = 0

    ratio[0] = vecteurs_Q[0] / matrice_coef_bloc_contrainte[0][colonne_valeur_entrante]

    for i in range(1, len(matrice_coef_bloc_contrainte)):
        if (
            vecteurs_Q[i] / matrice_coef_bloc_contrainte[i][colonne_valeur_entrante]
            < min_value
        ):
            min_value = (
                vecteurs_Q[i] / matrice_coef_bloc_contrainte[i][colonne_valeur_entrante]
            )
            index_min_value = i

        ratio[i] = (
            vecteurs_Q[i] / matrice_coef_bloc_contrainte[i][colonne_valeur_entrante]
        )

    coef_pivot[index_min_value] = min_value
    pivot = matrice_coef_bloc_contrainte[index_min_value][colonne_valeur_entrante]
    ligne_pivot = [x / pivot for x in matrice_coef_bloc_contrainte[index_min_value]]

    vecteurs_Q[index_min_value] = vecteurs_Q[index_min_value] / pivot

    matrice_coef_bloc_contrainte[index_min_value] = ligne_pivot

    # Calcul des lignes par rapport à la ligne pivot
    for row in range(len(matrice_coef_bloc_contrainte)):
        if row != index_min_value:
            coef_colonne_pivot = matrice_coef_bloc_contrainte[row][colonne_valeur_entrante]
            print("coef_colonne_pivot", coef_colonne_pivot)
            print("row de matrice_coef_bloc_contrainte", matrice_coef_bloc_contrainte[row])
            matrice_coef_bloc_contrainte[row] = [
                x - coef_colonne_pivot * pivot for x in matrice_coef_bloc_contrainte[row]
            ]
            vecteurs_Q[row] = vecteurs_Q[row] - (coef_colonne_pivot * pivot)

    print("ligne pivot", ligne_pivot)
    print("vecteurs_Q", vecteurs_Q)

    # division par le pivot pour la ligne pivot

    # Calcul de Z
    Z = 0
    for i in range(len(coef_pivot)):
        Z = Z + coef_pivot[i] * vecteurs_Q[i]

    # Calcul pour trouver le nouveau Valeur Dans la Base
    max_difference_Cj_Zj = max(difference_Cj_Zj)
    print("Max difference Cj - Zj: ", max_difference_Cj_Zj)
    index_max_difference_Cj_Zj = difference_Cj_Zj.index(max_difference_Cj_Zj)
    print("Index max difference Cj - Zj: ", index_max_difference_Cj_Zj)

    vdb[index_min_value] = "x" + str(index_max_difference_Cj_Zj + 1)

    # ============== Affichage du tableau Simplexe =================
    display_tab(
        coef,
        matrice_coef_bloc_contrainte,
        vecteurs_Q,
        coef_pivot,
        Z,
        Zj,
        difference_Cj_Zj,
        vdb,
        ratio,
    )
    # ===================== Fin de l'affichage =====================

    print("Pivot: ", pivot)
    print(
        "Position pivot: row:",
        index_min_value + 1,
        "column: x",
        colonne_valeur_entrante + 1,
    )


# Fonction pour afficher le tableau Simplexe
def display_tab(
    coef: List[float],
    matrice_sous_contrainte: List[List[float]],
    vecteurs_Q,
    coef_pivot,
    Z,
    Zj,
    difference_Cj_Zj,
    vdb,
    ratio,
):
    print("=================================================")
    print("Tableau Simplexe")
    print("Cj              |", coef, "| Ratio")
    print("-------------------------------------------------")
    if len(coef_pivot) == 2:
        print("CP   | vdb | Q  | x1 x2 x3 x4  |")
    else:
        print("CP   | vdb | Q  | x1 x2 x3 x4 x5 x6  |")
    print("-------------------------------------------------")
    for column in range(len(coef_pivot)):
        print(
            f"{coef_pivot[column]:.2f} | {vdb[column]}   | {vecteurs_Q[column]:.1f} | {matrice_sous_contrainte[column]} | {ratio[column]:.2f}"
        )
    print("Z =", Z, "   | Zj |", Zj)
    print("Cj - Zj         |", difference_Cj_Zj)


# =================== Testing ===================
# Coefficients de la fonction objectif
c = [4, 5, 6]

# Matrice des contraintes
A = [[2, 3, 1], [1, 4, 2], [3, 1, 2]]

# Côtés droits des inégalités
b = [10, 12, 14]

type_inegalite = ["<=", "<=", "<="]

# Appel de la fonction pour maximiser
print("3 Variables")
simplexe(c, A, b, type_inegalite, maximisation=True)

# Exemple avec 2 variables
c = [2, 3]
# Matrice des contraintes
A = [[1, 2], [2, 1]]

# Côtés droits des inégalités
b = [8, 10]

type_inegalite = ["<=", "<="]

# Appel de la fonction pour maximiser
# print("2 Variables")
# simplexe(c, A, b, type_inegalite, maximisation=True)
