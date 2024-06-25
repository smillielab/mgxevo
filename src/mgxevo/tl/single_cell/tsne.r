
run_largevis = function(data, k, save_knn=FALSE, save_weights=FALSE, dist.use='euclidean', verbose=TRUE, max_iter=1, perplexity=max(50, k/3)){
    require(igraph)
    require(largeVis)

    # Run largeVis on data with k nearest neighbors
    # save_knn saves the sparse knn graph
    # save_weights saves the sparse knn distance matrix

    # capitalize distance metric
    dist.use = paste(toupper(substr(dist.use, 1, 1)), substr(dist.use, 2, nchar(dist.use)), sep="")

    # run largevis
    res = largeVis(data, K=k, save_neighbors=save_knn, save_edges=save_weights, distance_method=dist.use, verbose=verbose, max_iter=max_iter, perplexity=perplexity)

    # construct output object
    out = list(coords=t(res$coords))
    if(save_knn == TRUE){
        # sparse KNN matrix -> KNN graph
        knn = res$knn
	knn = neighborsToVectors(t(knn))
	knn = graph_from_edgelist(as.matrix(data.frame(i=knn$i+1, j=knn$j+1)))
	out$knn = knn
    }
    if(save_weights == TRUE){
        # may need to convert to sparse matrix
	# e.g. sparseMatrix(i=, j=, x)
        out$weights = res$edges
    }
    return(out)
}

project_tsne = function(pca.rot.new, pca.rot.old, tsne.rot.old, perplexity, n.cores=1){
    source('~/code/single_cell/parallel.r')

    # Project new data onto old TSNE space
    sum_x_old = rowSums(pca.rot.old^2)

    # Get new TSNE rotations
    tsne.rot.new = run_parallel(foreach(i=rownames(pca.rot.new), .combine=cbind, .export=c('project_map', 'd2p_cell', 'KullbackFun', 'Hbeta')) %dopar% {
        project_map(pca.rot.new[i,], pca.rot.old, sum_x_old, tsne.rot.old)
    }, n.cores=n.cores)

    # Return [cells x TSNE axes] rotations
    colnames(tsne.rot.new) = rownames(pca.rot.new)
    return(as.data.frame(t(tsne.rot.new)))
}

project_map = function(z, x_old, sum_X_old, x_old_tsne, P_tol=5e-6, perplexity=30) {
  sum_z=sum(z^2)
  D_org=sum_z+(-2*as.matrix(x_old)%*%t(z)+sum_X_old)
  P=d2p_cell(D_org,perplexity)
  nn_points = which(P> P_tol);    #Use only a small subset of points to comupute embedding. This keeps all the points that are proximal to the new point
  X_nn_set = x_old[nn_points,]  #Original points
  y_nn_set = x_old_tsne[nn_points,];  #Computed embeddings
  P_nn_set = P[nn_points,]; #Probabilities
  y_new0 = (t(as.matrix(y_nn_set))%*%t(as.matrix(rbind(P_nn_set,P_nn_set))))[,1] #Initial guess of point as a weighted average
  sink("/dev/null")
  y_new=optim(y_new0,KullbackFun,gr=NULL,y_nn_set,P_nn_set,method = "Nelder-Mead")
  sink()
  return(y_new$par)
}

d2p_cell=function(D,u=15,tol=1e-4) {
  betamin=-Inf; betamax=Inf;
  tries = 0; tol=1e-4; beta=1
  beta.list=Hbeta(D,beta); h=beta.list[[1]]; thisP=beta.list[[2]]; flagP=beta.list[[3]]
  hdiff = h - log(u);
  while (abs(hdiff) > tol && tries < 50) {
    if (hdiff > 0) {
      betamin = beta;
      if (betamax==Inf) {
        beta = beta * 2;
      } else
      {
        beta = (beta + betamax) / 2;
      }
    } else {
      betamax = beta;
      if (betamin==-Inf) {
        beta = beta / 2;
      } else {
        beta = (beta + betamin) / 2;
      }
    }
    beta.list=Hbeta(D,beta); h=beta.list[[1]]; thisP=beta.list[[2]]; flagP=beta.list[[3]]
    hdiff = h - log(u);
    tries = tries + 1
  }
  # set the final row of p
  P=thisP
  #Check if there are at least 10 points that are highly similar to the projected point
  return(P)
}

KullbackFun=function(z,y,P) {

  #Computes the Kullback-Leibler divergence cost function for the embedding x in a tSNE map
  #%P = params{1};              %Transition probabilities in the original space. Nx1 vector
  #%y = params{2};              %tSNE embeddings of the training set. Nx2 vector

  print(z); print(dim(y))
  Cost0 = sum(P * log(P));          #Constant part of the cost function
  #Compute pairwise distances in embedding space
  sum_z=sum(z^2)
  sum_y=rowSums((y^2))
  D_yz=sum_z+(-2*as.matrix(y)%*%t(matrix(z,nrow=1))+sum_y)
  Q = 1/(1+D_yz);
  Q = Q/sum(Q)               #Transition probabilities in the embedded space
  Cost = Cost0 - sum(P*log(Q)); #% - 100 * sum(Q .* log(Q));
  return(Cost)
}

Hbeta=function(D, beta) {
    flagP=1;
    P = exp(-D * beta);
    sumP = sum(P);
    if (sumP<1e-8) { #In this case it means that no point is proximal.
      #P = ones(length(P),1) / length(P);
      P = matrix(rep(1, length(P))/length(P))
      sumP = sum(P);
      flagP=0;
    }
    H = log(sumP) + beta * sum(D * P) / sumP;
    P = P / sumP;
    return(list(H, P, flagP))
}
