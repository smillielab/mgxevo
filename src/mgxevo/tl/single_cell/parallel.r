require(doParallel)

bg.run = function(f, name=''){

   # Run function f in background
   p = mcparallel(f)

   # Check function f in background
   g = function(){
       # Wait for function f to finish
       q = mccollect(p, wait=TRUE)
       cat(paste0('\n\n', name, ' finished\n\n'))
       alarm()
       # Return function f output
       return(q[[1]])
   }
   mcparallel(g())
}

bg.get = function(p){
   mccollect(p)[[1]]
}

run_parallel = function(f, n.cores=1){

    # Execute parallel function on cluster
    # Example:
    # g = run_parallel(foreach(i=1:10, .combine=cbind), %dopar% rnorm(100, i), n.cores=10)

    if(n.cores > 1){

        # Register cluster
	cluster = makeCluster(n.cores, type='FORK', outfile='')
	registerDoParallel(cluster)

        # Run f in parallel
    	y = f

	# Register sequential
	stopCluster(cluster)
    	registerDoSEQ()

    } else {

	# Run f sequentially
        y = f
    }

    return(y)
}

p.saveRDS = function(object, file){

    # Create process to save object
    p = mcparallel(saveRDS(object, file))

    # Wait for saveRDS to complete
    f = function(){
        q = mccollect(p, wait=TRUE)
	print('Finished')
	alarm()
	return(q)
    }
    q = mcparallel(f())
}

