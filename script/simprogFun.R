
#Use a  paramfile===============================================================================================================
launchprogevol = function(myfile="param.txt") { #launch the program with a parameter file à mettre entre "" - ca fonctionne

  theparamsare=read.param(myfile)

  L=      gv("GENET_NBLOC", theparamsare)
  theta = gv("FITNESS_OPTIMUM", theparamsare)  ; if (length(theta)<L){ for (i in (length(theta)+1):L){theta[length(theta)+1]=theta[1]}} # on verifie que theta a le bon nombre de parametres, sinon on lui rajoute le 1er parametre le nombre de fois qu'il faut
  s1=     gv("FITNESS_STRENGTH", theparamsare) ; if (length(s1)<L){ for (i in (length(s1)+1):L){s1[length(s1)+1]=s1[1]}}
  s2=     gv("FITNESS_STABSTR", theparamsare)  ; if (length(s2)<L){ for (i in (length(s2)+1):L){s2[length(s2)+1]=s2[1]}}

  runevo=runevolution(N=gv("INIT_PSIZE", theparamsare),
                      gen=gv("SIMUL_GENER", theparamsare),
                      L=L,
                      mucis=gv("GENET_MUTRATES", theparamsare),
                      mutrans=gv("GENET_TRANSMUTRATES", theparamsare),
                      muteff=gv("GENET_MUTSD", theparamsare),
                      a=gv("INIT_BASAL", theparamsare),
                      theta=theta,
                      s1=s1,
                      s2=s2,
                      pas=gv("SIMUL_OUTPUT", theparamsare),
                      steps=gv("DEV_TIMESTEPS", theparamsare),
                      measure=gv("DEV_CALCSTEPS", theparamsare),
                      diag=gv("INIT_CONDIAG", theparamsare),
                      clon=gv("INIT_CLONAL", theparamsare))
  return(runevo)
} #end launchevol2

gv=function(input, param2){
	# gv: get value
	# First parameter:

  stopifnot(is.character(input), length(input)==1)
  stopifnot(is.list(param2))
  if (! input %in% names(param2)) stop("Tag ", input, " not found in the parameter file.")

  return(param2[[input]])
} #end gv

gv.test=function(input,myfile="param.txt"){dat=read.param(myfile);value = as.numeric(unlist(strsplit(dat[[input]],",")));return(value)}

read.param = function(paramfile) { #from Estelle: turn un fichier genre table into a list
	filterpar = function(line) {
		line = line[-1] # Remove the name tag

		# Two possibilities:
		#    the line contains only numbers -> returns a numeric vector
		#    the line contains at least a non-number -> returns a character vector
		ans = suppressWarnings(as.numeric(line))
		if(any(is.na(ans))) ans = line
		return(ans)
	}

	stopifnot(is.character(paramfile), file.exists(paramfile))
	ss  = scan(file=paramfile, what=character(), sep="\n", quiet=TRUE)
	ss  = strsplit(ss, split="\\s+")
	nam = sapply(ss, "[", 1)
	ans = lapply(ss, filterpar)
	names(ans) = nam
	return(ans)
} #end readparam

#===================================================================================================================
#===================================================================================================================

runevolution = function(N=10, gen=10, L=3, mucis=0.5, mutrans=0.5, a=0.2, pas=1, steps, measure, diag, muteff, clon, ...){ ##N=individus gen=generations L=genes ...qui permettent d'acceder aux parametres des functions downstream

  if (clon=="clonal"){    #~~mettre en init_clonal yesss ##hard
        pop = development(populclonale(N, L,diag, ...), N, a, L, steps, measure,...)
  }else{pop = development(populrandom(N, L, diag, ...), N, a, L, steps, measure,...)}

  #create and name columns of the output dataframe
  output = data.frame(matrix(ncol = 1+L*6+2+L*L*2, nrow = 0))
  names=c("Gen")
  for (i in c("MPhen","VPhen","MUnstab","VUnstab", "MTrans", "VTrans")){
    for (j in 1:L){names[length(names)+1]=paste(i, j, sep="")}}
  names[length(names)+1] = "MFit"
  names[length(names)+1] = "VFit"
  for (i in c("MeanAll","VarAll")){
    for (j in 1:square(L)){names[length(names)+1]=paste(i, j, sep="")}}
  colnames(output) = names

  count1=0
  for (gg in 1:gen){
    #SORTIE DES DATAS POUR ANALYSES
    if(gg==1 || gg == gen || gg %% pas==0) { #definition du pas pour sortir les datas, pour avoir la premiere, la derniere et sorties intermediaires
      output[nrow(output)+1,] = NA #on rajoute une ligne vide
      output[nrow(output), "Gen"] = gg
      output[nrow(output), "MFit"] = mean(sapply(1:N, function(i) pop[[i]]$fitness))
      output[nrow(output), "VFit"] = var(sapply(1:N, function(i) pop[[i]]$fitness))
      for (i in 1:L){
        output[nrow(output), paste("MPhen",   i, sep="")] = mean(sapply(1:N, function(j) pop[[j]]$mean[[i]]))
        output[nrow(output), paste("VPhen",   i, sep="")] =  var(sapply(1:N, function(j) pop[[j]]$mean[[i]]))
        output[nrow(output), paste("MUnstab", i, sep="")] = mean(sapply(1:N, function(j) pop[[j]]$var[[i]]))
        output[nrow(output), paste("VUnstab", i, sep="")] =  var(sapply(1:N, function(j) pop[[j]]$var[[i]]))
        output[nrow(output), paste("MTrans",   i, sep="")] = mean(sapply(1:N, function(j) pop[[j]]$ind[[i,L+1]]))
        output[nrow(output), paste("VTrans",   i, sep="")] =  var(sapply(1:N, function(j) pop[[j]]$ind[[i,L+1]]))
      }
      for (i in 1:square(L)){
        output[nrow(output), paste("MeanAll", i, sep="")] = mean(sapply(1:N, function(j) pop[[j]]$ind[i]))
        output[nrow(output), paste("VarAll",  i, sep="")] =  var(sapply(1:N, function(j) pop[[j]]$ind[i]))
      }
    } #end of sortie

    pop = development(newpopul(pop, N, L, mucis, mutrans, diag, muteff), N, a, L, steps, measure, ...)

    #cat("\014"); cat("\n"); cat("",paste(Sys.time()),"\n")
    #cat(" RUNNING Generation ", str_pad(gg,3,pad=" ")," - ", str_pad(round(gg/gen*100),3,pad=" "),"%  ", dna(floor(gg/gen*100/2)),rep("-",(50-floor(gg/gen*100/2))),"\n", sep="")#, dna(2*gg%%25)
    count2=round(gg/gen*100/2)
    ##if (count1<count2){count1=count2; cat(dna(1))}
  } #end of generation

  # cat("\n")
  #cat(" SUCCESSFULLY COMPLETED, With a ")
  #return(pop)
  return(output)
}

#POP================================================================================================================
#===================================================================================================================

    #######    ##     ##    ####      ###   #F#####   #F#########  ###    #######     ####      ###    #######
  #FF#######  ###     ###  ######     ###  #FF######  F##########  ###  ###########  ######     ###  ###########
  #FF    ###  ###     ###  #######    ###  #FF    ##      #FF      ###  ###     ###  #######    ###  ###     ###
  #FF         ###     ###  ###   ##   ###  #FF            #FF      ###  ###     ###  ###   ##   ###  ###
  #FF######   ###     ###  ###    ### ###  #FF            #FF      ###  ###     ###  ###    ### ###  #####
  #FF######   ###     ###  ###     ######  #FF            #FF      ###  ###     ###  ###     ######    #######
  #FF         ###     ###  ###      #####  #FF            #FF      ###  ###     ###  ###      #####        #####
  #FF         ###     ###  ###       ####  #FF            #FF      ###  ###     ###  ###       ####          ###
  #FF         ###     ###  ###        ###  #FF     ##     #FF      ###  ###     ###  ###        ###  ###     ###
  #FF         ###########  ###        ###  #FF#######     #FF      ###  ###########  ###        ###  ###########
  #FF           #######    ###        ###   #F######      #FF      ###    #######    ###        ###    #######

#POP================================================================================================================
#===================================================================================================================

indiv = function(L=5, diag=1, ...){
  ind=list()                   #INIT_ALLELES
  ind$mom = matrix(rnorm(L*L,0,0.1), nrow=L) #matrix(rep(0,square(L)), nrow=L) ind$mom =  #haplo mother #Commence avec des matrix vides pour checker si ca marche indeed
  ind$dad = matrix(rnorm(L*L,0,0.1), nrow=L) #matrix(rep(0,square(L)), nrow=L) ind$dad =  #haplo father
  if (diag==0){for (i in 1:L){ind$mom[i,i]=0 ; ind$dad[i,i]=0}}
  trans=rep(1,L)
  ind$mom = cbind(ind$mom,trans)
  ind$dad = cbind(ind$dad,trans)
  ind$ind = (ind$mom+ind$dad)/2
  return(ind)
}

indivhomoz = function(L=5, diag=1, ...){
  ind=list()                   #INIT_ALLELES
  ind$mom = matrix(rnorm(L*L,0,0.1), nrow=L) #matrix(rep(0,square(L)), nrow=L) ind$mom =  #haplo mother #Commence avec des matrix vides pour checker si ca marche indeed
  ind$dad = ind$mom
  if (diag==0){for (i in 1:L){ind$mom[i,i]=0 ; ind$dad[i,i]=0}}
  trans=rep(1,L)
  ind$mom = cbind(ind$mom,trans)
  ind$dad = cbind(ind$dad,trans)
  ind$ind = (ind$mom+ind$dad)/2
  return(ind)
}

populrandom = function(N=10, L=5, diag=1, ...){ #ind, #genes
  pop=list()
  #for (i in 1:N) {pop[[i]]=indiv(L, diag, ...)}
  pop = lapply(1:N, function(i) pop[[i]]=indiv(L, diag, ...))
  return(pop)
}

populclonale = function(N=10, L=5, diag=1, ...){
  pop = list()
  Lelu = indivhomoz(L, diag,...) #car cest lui l'elu
  #for (i in 1:N) {pop[[i]]=Lelu} #OPTIMISATION:apply
  pop = lapply(1:N, function(i) pop[[i]]=Lelu)
  return(pop)
}

populdevrandom = function(){
  test=development(populrandom(N=10, L=5, diag=0))
  #for (i in 1:10){test[[i]]$fitness=runif(1)}
  return(test)
}

#DEV================================================================================================================
#===================================================================================================================
development = function(pop, N=10, a=0.2, L=5, steps=40, measure=10, ...){
  pop = lapply(pop, function(i){ #~mclapply
    dev = model.M2(i$ind, a, steps, measure)
    i$mean = dev$mean
    i$var  = dev$var
    i = addfitness(i, L, ...)
    })#, mc.cores=3)
  return(pop)
}

sigma.M2p = function(x, aam1, l1am1) { 1. / (1. + exp((-x/aam1)+l1am1)) }
model.M2  = function(W, a=0.2, steps=40, measure=10, init= rep(a,nrow(W)), S0=init, full=FALSE, varFUN=function(x) mean((x-mean(x))^2)) {
#note# steps:pas de temps du developement  #measure:on fait les mesure sur les x derniers steps du development
#note# init: runif(nrow(W), min=0, max=1), ##ici toutes les valuers initiales sont aleatoires : entre individus et entre generaztions...attention
#note# S0:#minimum,maximum,median,random_binary,random,basal
  aam1 = a*(1-a) ; l1am1 = log(1/a-1)
  sto = matrix(NA, nrow=length(S0), ncol=steps+1)
  sto[,1] = S0
  L=nrow(W)
  transeffect=W[,(L+1)] ; W=W[,1:L]  #on separe les arguments qui vont passer dans l'equation
  W2 = t(t(W)*transeffect)
  for (i in 1:steps) {
    S0 = sigma.M2p(W2 %*% S0, aam1, l1am1) 	#attention, et pas S0 %*% W
    sto[,i+1] = S0
  }
  ans = list()
  ans$mean = apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
  ans$var  = apply(sto[,(steps+1-measure):(steps+1)], 1, varFUN)
  if (full) ans$full = sto
  return(ans)
}

square=function(n){q=n*n;return(q)}

addfitness = function(ind, L=5, theta=rep(1,L), s1=rep(10,L), s2=rep(4600,L), ...) {#optimum pour les 3 genes a la fois et #force de selection
  fitmean = 0 ; fitvar = 0
  for (i in 1:L) {
    fitmean = fitmean + (-s1[i]*(square(ind$mean[i]-theta[i])))
    fitvar  = fitvar  + (-s2[i]*(ind$var[i]))
    }
  ind$fitness = exp(fitvar)*exp(fitmean)
  return(ind)
}

#NEWPOP=============================================================================================================
#===================================================================================================================

newpopul = function(pop, N=10, ...){
  newpop = list()
  newpop = lapply(1:N, function(i) newpop[[i]]=reproduction(pop, N, ...))
  #for (i in 1:N) {newpop[[i]]=reproduction(pop, N, ...)}
  return(newpop)
}

reproduction = function(pop, N, ...){#create new individu from pop
  ind=list()
  #hermaph=TRUE ; while (hermaph==TRUE) { #Control autofecondation
    mom = selectforfitness(pop, N)
    dad = selectforfitness(pop, N)
  #  if (mom != dad){hermaph=FALSE}
  #}
  ind$mom = gametogenesis(pop[[mom]], ...)
  ind$dad = gametogenesis(pop[[dad]], ...)
  ind$ind = (ind$mom+ind$dad)/2
  return(ind)
}

selectforfitness.old= function(pop, N=10){ #select one individual proportionally to fitness
  fit=matrix(NA, N, 2)
  for (i in 1:N){
    fit[i,1]=i
    fit[i,2]=pop[[i]]$fitness
  }
  fit[,2]=cumsum(fit[,2])
  fit[,2]= fit[,2]/max(fit[,2])
  r=runif(1)
  fit[,2]=abs(fit[,2]-r)
  fit=fit[order(fit[,2]),] #he oui on choisit celle avec le plus petit ecart
  return(fit[1,1])
} # PROBLEM IN THIS FUNCTION : THE PROBABILITIES FOR EACH INDIVIDUAL ACCORDING TO ITS FITNESS DONT STRICTLY REFLECT THE EXPECTED ONES

selectforfitness = function(pop, N=10) {
	fitnesses <- sapply(pop, "[", "fitness")
	return(sample.int(length(pop), 1, prob=fitnesses))
}

gametogenesis = function(ind, L=5, mucis=0.002, mutrans=0.001, diag=1, muteff, ...){ #ind1=test[[10]] #from ind, create gamete: a matrix
  gamete=matrix(NA, L,(L+1))
  mutdad=ind$dad ; mutmom=ind$mom
  mutdad=cismutation(mutdad, L, mucis, diag, muteff)     ; mutmom=cismutation(mutmom, L, mucis, diag, muteff)
  mutdad=transmutation(mutdad, L, mutrans, diag, muteff) ; mutmom=transmutation(mutmom, L, mutrans, diag, muteff)
  for (i in 1:L){
    r=runif(1)
    if (r<0.5){gamete[i,]=mutdad[i,]}else{gamete[i,]=mutmom[i,]}   ##REC TAUX DE RECOMBINAISON GENET_RECRATES
  } #Attention#~ taux de recombinaison marche pas pour different de 0.5 car on sattend a ce que ci recombinaison on reste ensuite sur mom ou dad, sur lequel on a changE
  return(gamete)
}

cismutation = function(W, L=5, mucis=0.001, diag=1, muteff){
  transeffect=W[,L+1] ; W=W[,1:L]
  r=rpois(1,mucis)
  if (r>0){ # on tire un nombre de mutation sur une loi de poisson - on utilise si egal a 1 ou si superieur a 1???
    for (i in 1:r){ #car il peut tirer r=2 donc 2 mutations sur l'individu
      gene=sample(c(1:L),1) ; bindsite=sample(c(1:L),1)#tire quel gene va etre mutE, quelle ligne
      if (diag==0){while (gene==bindsite){ bindsite=sample(c(1:L),1) }} #on garde les diag à zero
      combien=rnorm(1,0,muteff) #on tire de combien on va muter
      W[bindsite,gene]=W[bindsite,gene]+combien
    }
  }
  W=cbind(W,transeffect)
  return(W)
}

transmutation = function(W, L=5, mutrans=0.001, diag=1, muteff){ #for the real trans effect on the vector
  transeffect=W[,L+1] ; W=W[,1:L]
  r=rpois(1,mutrans)
  if (r>0){
    for (i in 1:r){
      gene=sample(c(1:L),1)
      combien=rnorm(1,0,muteff)
      transeffect[gene]=transeffect[gene]+combien
      if (transeffect[gene]>2){transeffect[gene]=2}
      if (transeffect[gene]<0){transeffect[gene]=0}
    }
  }
  W=cbind(W,transeffect)
  return(W)
}

#GRAPHS#######################################################################################################
#Traitement des datas from dataframe#######################################################################################################

#UN PEU PLUS VIF:
bgcol=function(){bg='#e7ebf4'; return(bg)} #899caa
mycolors=function(){mycolors=c('#b53122','#892daf','#2c85c2','#25bd66','#e09b2d','#dd5c25','#0b0b0a','#cd30cc','#29cb96') ; return(mycolors)}
varcol=function(){varcol='#dbdbed'; return(varcol)} #8092a6

##DOESNT HURT MY EYES:
#bgcol=function(){bg='#899caa'; return(bg)}
#mycolors=function(){mycolors=c('#78281f','#4a235a','#1b4f72','#186a3b','#b9770e','#b9410e','#0b0b0a','#9e3c9d','#3c9e7e') ; return(mycolors)}
#varcol=function(){varcol='#8092a6'; return(varcol)}


graphfitexp = function(output=test){ #output est en fait l'input ici, cad le dataframe
  par(mfrow=c(2,1), bg=bgcol())
  graphfit(output)
  graphexp(output)
}

graphfit = function(output=test){
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  par(bg=bgcol())
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(0,1), main="Fitness of the population", xlab="Generations", ylab="Fitness of pop")
  polygon(c(output$Gen, rev(output$Gen)), c(output$MFit+output$VFit, rev(output$MFit-output$VFit)), col=varcol(), border=NA)
  lines(output$Gen, output$MFit, type = 'l', col='#001f3f')
}

graphexp=function(output=test){
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  par(bg=bgcol())
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(0,1), main="Gene expressions - Phenotype moyen", xlab="Generations", ylab="Expression")
  mycolors=mycolors() ; if (length(mycolors)<L) {for (i in (length(mycolors)+1):L){mycolors[i]=mycolors[1]}} #si ya pas assez de couleurs, ce sera affiché quand meme avec la premiere couleur
  for (i in 1:L){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MPhen", i, sep="")]]+output[[paste("VPhen", i, sep="")]], rev(output[[paste("MPhen", i, sep="")]]-output[[paste("VPhen", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:L){ lines(output$Gen, output[[paste("MPhen", i, sep="")]], type = 'l', col=mycolors[i])}
  legend("topleft", lty=1, col=mycolors, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphmtrxvalues = function(output=test){ #GRAPH les valeurs de la matrice
  par(bg=bgcol())
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  tt=output[,grep(x = colnames(output), pattern = "MeanAll.*")] #extracts values of the MeanAll
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(min(tt)-0.1,max(tt)+0.1), main="Mean Matrix values (alleles) in population", xlab="Generations", ylab="Strength of interaction")
  mycolors=mycolors() ; newcolors=c();
  #colors par colonnes (par effet DU gene sur les autres):
  for (j in 1:L){for (i in 1:L){newcolors[length(newcolors)+1]=mycolors[j]}} ; mtext("Effect OF the colored gene on others", 3, line=.3)
  #colors par lignes (par effet SUR gene des autres):
  #for (j in 1:L){for (i in 1:L){newcolors[length(newcolors)+1]=mycolors[i]}} ; mtext("Effect ON the colored gene by others", 3, line=.3)
  for (i in 1:(L*L)){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MeanAll", i, sep="")]]+output[[paste("VarAll", i, sep="")]], rev(output[[paste("MeanAll", i, sep="")]]-output[[paste("VarAll", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:(L*L)){ lines(output$Gen, output[[paste("MeanAll", i, sep="")]], type = 'l', col=newcolors[i])}
  legend("topleft", lty=1, col=mycolors, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphtransvector = function(output=test){
  par(bg=bgcol())
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  tt=output[,grep(x = colnames(output), pattern = "MTrans.*")]
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(0,2), main="Trans values in pop", xlab="Generations", ylab="Strength of interaction / activation force") #ylim=c(min(tt)-0.1,max(tt)+0.1)
  mycolors=mycolors() ; if (length(mycolors)<L) {for (i in (length(mycolors)+1):L){mycolors[i]=mycolors[1]}}
  for (i in 1:L){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MTrans", i, sep="")]]+output[[paste("VTrans", i, sep="")]], rev(output[[paste("MTrans", i, sep="")]]-output[[paste("VTrans", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:L){ lines(output$Gen, output[[paste("MTrans", i, sep="")]], type = 'l', col=mycolors[i])}
  legend("topleft", lty=1, col=mycolors, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphall=function(output=test){
  par(mfrow=c(2,2), bg=bgcol())
  graphfit(output)
  graphmtrxvalues(output)
  graphexp(output)
  graphtransvector(output)
}

matrixFinale.Matrx = function(fichierdesortie=test){ #extraire la matrice moyenne de la derniere generation du dataframe de sortie
  L=length(grep(x = colnames(fichierdesortie), pattern = "MPhen.*")) #guess L from number of columns of the dataframe
  taille=nrow(fichierdesortie) #test$gen[nrow(test)]
  lastW=c()
  for (i in 1:(L*L)){ lastW[length(lastW)+1] = fichierdesortie[[paste("MeanAll", i, sep="")]][taille] }
  for (i in 1:L)    { lastW[length(lastW)+1] = fichierdesortie[[paste("MTrans",  i, sep="")]][taille]  }
  lastW=matrix(lastW, nrow=L)
  transeffect=lastW[,(L+1)] ; lastW=lastW[,1:L]
  lastW=t(t(lastW)*transeffect)
  return(lastW)
}

matrixFinale.Genotype = function(fichierdesortie=test){ #extraire la matrice moyenne de la derniere generation du dataframe de sortie
  L=length(grep(x = colnames(fichierdesortie), pattern = "MPhen.*"))#guess L from number of columns of the dataframe
  taille=nrow(fichierdesortie) #test$gen[nrow(test)]
  lastW=c()
  for (i in 1:(L*L)){ lastW[length(lastW)+1] = fichierdesortie[[paste("MeanAll", i, sep="")]][taille] }
  for (i in 1:L)    { lastW[length(lastW)+1] = fichierdesortie[[paste("MTrans",  i, sep="")]][taille]  }
  lastW=matrix(lastW, nrow=L)
  return(lastW)
}

#==============================================================================================
#==============================================================================================
#==============================================================================================
#==============================================================================================

dna=function(n=2500){ #just for fun
  paste(sample(c("A", "T", "G", "C"), n, replace=TRUE), collapse = '')
}
