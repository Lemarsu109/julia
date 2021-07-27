########## Bibliothèques ##########
using BioStructures
using Bio3DView
using Blink
using BioAlignments
using BioStructures: read as reader
using Plots: plot, plot!
using LinearAlgebra: svd, det


########## Fonctions de lecture ##########
# Lecture du PDB pour retourner une matrice de coordonnées (pour une protéine)
function read_pdb(path::String)
    struc = reader(path, BioStructures.PDB)
    # calphaselector permet de ne sélectionner que les carbones alpha
    L = collectatoms(struc, calphaselector)
    M = Array{Float64,2}(undef, length(L), 3) # Introduit une matrice vide
    for i = 1:length(L)
        M[i,:] = BioStructures.coords(L[i])
    end
    return M
end

# Lecture du PDB pour retourner une matrice de coordonnées (pour un modèle)
function read_pdb(path::String, model::Int64)
    struc = reader(path, BioStructures.PDB)
    # calphaselector permet de ne sélectionner que les carbones alpha
    L = collectatoms(struc[model]['A'], calphaselector)
    M = Array{Float64,2}(undef, length(L), 3) # Introduit une matrice vide
    for i = 1:length(L)
        M[i,:] = BioStructures.coords(L[i])
    end
    return M
end

# Lecture du PDB pour retourner une matrice de coordonnées (pour une chaîne)
function read_pdb(path::String, model::Int64, chain::Char)
    struc = reader(path, BioStructures.PDB)
    # calphaselector permet de ne sélectionner que les carbones alpha
    L = collectatoms(struc[model][chain], calphaselector)
    M = Array{Float64,2}(undef, length(L), 3) # Introduit une matrice vide
    for i = 1:length(L)
        M[i,:] = BioStructures.coords(L[i])
    end
    return M
end


########## Fonctions de visualisation ##########
# Permet de visualiser une protéine avec Bio3DView
function visual(path::String)
    struc = read(path, PDB)
    viewstruc(struc[1])
end


########## Fonctions d'alignements ##########
function align(struc1, struc2)
    seq1, seq2 = LongAminoAcidSeq.([struc1, struc2], standardselector, gaps=false)
    scoremodel = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)
    pairalign(GlobalAlignment(), seq1, seq2, scoremodel)
end


########## Fonctions de RMSD ##########
# Calcul le centroide d'une matrice
function centroid(mat)
    cx, cy, cz = 0, 0, 0
    for i = 1:size(mat, 1)
        cx += mat[i,1]
        cy += mat[i,2]
        cz += mat[i,3]
    end
    cx /= size(mat, 1)
    cy /= size(mat, 1)
    cz /= size(mat, 1)
    return cx, cy, cz
end

# Centre la matrice à l'origine du repère
function center_th!(mat2)
    cx2, cy2, cz2 = centroid(mat2)
    for i = 1:size(mat2, 1)
        mat2[i,1] -= cx2
        mat2[i,2] -= cy2
        mat2[i,3] -= cz2
    end
    return mat2
end

# Tourne la matrice à l'aide de l'algorithme de Kabsch
function rot_kabsch(mat1, mat2)
    matc = transpose(mat1) * mat2 # Matrice de covariance
    F = svd(matc) # Décomposition en valeurs singulières
    d = sign(det(F.V * transpose(F.U)))
    return F.V * [1 0 0; 0 1 0; 0 0 d] * transpose(F.U) # Matrice de rotation
end

# Calcul du RMSD à la main
function rmsd_alt(mat1, mat2)
    res = 0
    for i = 1:size(mat1, 1)
        res += (mat1[i,1] - mat2[i,1])^2
        res += (mat1[i,2] - mat2[i,2])^2
        res += (mat1[i,3] - mat2[i,3])^2
    end
    return sqrt(res/size(mat1,1))
end

# Fonction principale du RMSD à la main (2 protéines)
function rmsd_main(path1::String, path2::String)
    mat1 = read_pdb(path1)
    mat2 = read_pdb(path2)
    # Permet de mettre les matrices à la même taille
    #mat1 = mat1[1:size(mat2, 1),:]

    mat2 = center_th!(mat2)
    mat1 = center_th!(mat1)
    matr = rot_kabsch(mat1, mat2)
    mat2k = mat2 * matr

    println("RMSD fait main : ", rmsd_alt(mat1, mat2k))
    return mat1, mat2k
end

# Fonction principale du RMSD à la main (2 modèles d'une protéine)
function rmsd_main(path::String, model1::Int64, model2::Int64)
    mat1 = read_pdb(path, model1)
    mat2 = read_pdb(path, model2)
    # Permet de mettre les matrices à la même taille
    #mat1 = mat1[1:size(mat2, 1),:]

    mat2 = center_th!(mat2)
    mat1 = center_th!(mat1)
    matr = rot_kabsch(mat1, mat2)
    mat2k = mat2 * matr

    println("RMSD fait main : ", rmsd_alt(mat1, mat2k))
    return mat1, mat2k
end

# Fonction principale du RMSD à la main (2 chaînes de 2 modèles d'une protéine)
function rmsd_main(path::String, model1::Int64, model2::Int64, chain1::Char, chain2::Char)
    mat1 = read_pdb(path, model1, chain1)
    mat2 = read_pdb(path, model2, chain2)
    # Permet de mettre les matrices à la même taille
    #mat1 = mat1[1:size(mat2, 1),:]

    mat2 = center_th!(mat2)
    mat1 = center_th!(mat1)
    matr = rot_kabsch(mat1, mat2)
    mat2k = mat2 * matr

    println("RMSD fait main : ", rmsd_alt(mat1, mat2k))
    return mat1, mat2k
end

# Calcul du RMSD par BioStructures
function rmsd_julia(struc1, struc2)
    at1 = collectatoms(struc1)[calphaselector.(collectatoms(struc1))]
    at2 = collectatoms(struc2)[calphaselector.(collectatoms(struc2))]
    println("RMSD fait par BioStructures : ", BioStructures.rmsd(at1, at2))
end


########## Fonctions de graphes ##########
# Affiche le graphe des deux protéines en trois dimensions
function ploter(mat1, mat2, label1::String, label2::String)
    plot(mat1[:,1], mat1[:,2], mat1[:,3], label = label1, alpha = 0.5)
    plot!(mat2[:,1], mat2[:,2], mat2[:,3], label = label2, alpha = 0.5)
end

# Affiche la carte de contact d'une structure selon une distance
function mapcontact(struc, dist)
    contacts = ContactMap(collectatoms(struc, calphaselector), dist)
    plot(contacts)
end


########## Fonction principale (à décommenter) ##########
function main()
    # On donne les chemins des PDB et on les place dans des variables
    path = "/Users/aidanbonnefond/Downloads/pembrolizumab_models.pdb"
    struc = reader(path, BioStructures.PDB)
    p = "/Users/aidanbonnefond/Downloads/1SSU.pdb"
    s = reader(p, BioStructures.PDB)

    # On visualise la molécule
    #visual(path)
    #visual(p)

    # On aligne les séquences des protéines
    #print(align(struc[1], struc[2]))
    #print(align(s[1], s[2]))

    # On calcule le RMSD et on garde les matrices de coordonnées
    #rmsd_julia(struc[1]['A'], struc[2]['A'])
    #mat1, mat2 = rmsd_main(path, 1, 2, 'A', 'A')
    #rmsd_julia(s[2], s[3])
    #m1, m2 = rmsd_main(p, 2, 3)

    # On affiche les matrices comme alignement structural
    #ploter(mat1, mat2, "model 1", "model 2")
    #ploter(m1, m2, "model 1", "model 2")
    
    # On affiche les cartes de contacts
    #mapcontact(struc[1], 10.0)
    #mapcontact(s[1], 7.0)
    
end

main()