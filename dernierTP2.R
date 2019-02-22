ex4 = read.table("exo4.txt", sep = '\t', header = TRUE)
ex4

# Pop souris 
# Echantillon 1 N = 20
# Unite stats : un individu
# Exp aléa = traitement
# VA. P : "Poids"      }va quantitative
#     T : [H]          }
#
#Pop. E(X) = mu

#H0:rho=0 pas de lien   
#H1:rho/=0 lien

#-> Pearson -> lien linéaire      paramétrique
#-> Spearmann -> lien monotone    non paramétrique

Poids = ex4$Poids
Taux = ex4$Taux
cor.test(Poids,Taux)

cor(x = Poids, y = Taux)


cor.test(x = Poids,y = Taux,
         method = "spearman")


#R = Cov(X,Y)/(Sx * Sy) = 0.74

#t alpha=5%, (n-2)ddl =  2,101


data.frame(Taux,
           RgT = rank(Taux),
           Poids,
           RgP = rank(Poids),
           Diff = rank(Taux)-rank(Poids),
           Diffaucarre = (rank(Taux)- rank(Poids))^2)

model = lm(Taux~Poids)
model
plot(model)



#Spearman 
#Rcal= 1- 6somme (di^2/n(n-1)) = 0.7083

#Val théo = R (1-alpha/2), 20 ddl = 0.4451

#Rcal > Rthéo  donc H0 rejeté 
#il y a donc un lien




#      APRES
    



#y = ax + b 
#y = 0.45 * (-17.62)

#a^ = (Sxy)/(Sx²) = 14.03/30.51 = 0.46

#b^ = my - a*mx = 13.55 - 0.45 * 67.75
