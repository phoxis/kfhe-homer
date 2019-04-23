require (caret);
require (PerfMeas);
require (utiml);

###########################################
####             Training             #####
####################################################################################################
#### Function  : kfhe_ml_train                                                                  ####
#### Arguments :                                                                                #### 
####             X                : A matrix representing the input dataset                     ####
####             Y                : A matrix representing the target classes                    ####
####             max_iter         : Maximum ensemble iterations                                 ####
####             blend_with_class : If TRUE,  then the measurement updates will be done after   ####
####                                          converting to classes.                            ####
####                                If FALSE, then the measurement updates will be done on the  #### 
####                                          predicted scores                                  ####
####             early_stop       : If TRUE,  stop training after the Kalman gain is 1          ####
####                                If FALSE, continue training until 'max_iter' ensemble       ####
####                                          iterations                                        ####
####             reset_dist       : If TRUE,  then reset the training weights to uniform if     ####
####                                          the measurement returns error more                ####
####                                          than '1 / (1 - nclass)'                           ####
####                                If FALSE, do not reset the training weights                 ####
####             feedback         : Print the states of the parameters while training           ####
####                                                                                            ####
#### Return    : An object of type kfhe_m. Consists of kf-m. Other components are for           ####
####             debugging and analysis purpose                                                 ####
####################################################################################################
kfhe_ml_train <- function (X, Y, max_iter, method = c ("homer", "SVM"), blend_with_class = TRUE, early_stop = TRUE,  reset_dist = TRUE, scale_flag = FALSE, save_in_disk = FALSE, cores = 1, feedback = TRUE)
{
  if (feedback == TRUE) { cat ("submethod = (", method[1], ", ", method[2], ") \n\n"); }
  print_prec <- 8;  # Number of decimal places to round when printing feedback
  
  # Kalman filter variables for the model Kalman Filter, kf-m
  kf_m             <- list (); 
  kf_m$P           <- c ();    # Posterior variance/error
  kf_m$V           <- c ();    # Measurement variance/error
  kf_m$K           <- c ();    # Kalman Gain
  kf_m$h_list      <- list (); # List of learnt models
  kf_m$init_h      <- c ()     # The initial model.
  kf_m$train_blend <- c ();    # Internal state
    
  # Kalman filter variables for the kf-w
  kf_d        <- list ();
  kf_d$P_dist <- c ();    # Posterior variance/error
  kf_d$V_dist <- c ();    # Measurement variance/error
  kf_d$K_dist <- list (); # Kalman Gain
  kf_d$D      <- c ();    # Distribution
  
  # Debugging information
  debug_info                 <- list ();
  debug_info$train_accuracy  <- c ();    # Accuracy on the train blend error
  debug_info$train_f         <- c ();    # F-Score  on the train blend error
  debug_info$state_err       <- list (); # Error function evaluation of the train blend with respect to ground truth
  debug_info$blend_err_arr   <- c ();    # Uniformly weighted training blend error
  debug_info$D               <- list (); # Store distribution of all iterations
  debug_info$reset_pt        <- rep (FALSE, max_iter); # Trace in which iteration a distribution reset occurred
  debug_info$blend_err_delta <- c();
  
  
  cluster_types <- c ("balanced", "clustering", "random");
  kernel_types  <- c ("radial", "linear");# "polynomial");

  cluster_numbers <- c ();
  cluster_numbers <- c (2: sqrt (ncol (Y))); # WARNING
  if (feedback == TRUE) { cat ("cluster_number_ranges = ", cluster_numbers, "\n"); }

  
  # Final model to pack, stamp and return
  retmodel <- list ();
  
  while (1)
  {
    this_homer_cluster_method <- "balanced";
    kf_m$init_h      <- try (homer (mldr_from_dataframe (cbind (X, Y),  labelIndices=(ncol(X)+1):(ncol(X)+ncol(Y))), cores = cores, base.algorithm = method[2], method = this_homer_cluster_method, clusters = cluster_numbers[sample(length(cluster_numbers))[1]], scale = scale_flag), TRUE);
    
    if (class (kf_m$init_h) != "try-error")
    {
      break;
    }
    this_homer_cluster_method <- "random";
    cat ("Handling error\n");
    print (kf_m$init_h);
  }
  kf_m$train_blend <- predict (kf_m$init_h, X);
 
 
  unwt_comp_err    <- err_fun (kf_m$train_blend, Y, NULL); # Find the per datapoint error. This is a vector.
  uniwt_comp_err   <- unwt_comp_err * (1/nrow (X));        # Get a uniformly weighted version. This is a weighted vector.
  
  # Initialise the state variance for kf-m
  kf_m$P[1] <- 1.0; # No confidence
  
  # Initialise the state variance for kf-w
  # Optionally we can consider this as a vector and consider each component of
  # the distribution being managed by the KF
  kf_d$D              <- rep (1 / nrow (X), nrow (X));
  kf_d$P_dist         <- 1.0;
  
  debug_info$D[[1]]   <- kf_d$D; #initial training distribution
  
  this_blend_err      <- sum (err_fun(kf_m$train_blend, Y, NULL) * (1/nrow (X)));
  
  unique_number <- round (rnorm (1) * 10000);
  
  ###################################
  #### Start ensemble iterations ####
  ###################################
  for (t in 1:max_iter)
  {
    # Print feedback
    if (feedback == TRUE) { cat ("Iter = ", formatC (t, 4), ", "); }
    ############################################
    ############################################
    
    
    ############################
    #### Time update kf-m   ####
    ############################
    proc_noise <- 0;
    kf_m$train_blend <- kf_m$train_blend;
    kf_m$P[t]        <- kf_m$P[t] + proc_noise;
    
    #################################
    #### Measurement update kf-m ####
    #################################
    
    # Measurement process
    # Retry loop, if the error is more than a threshold then recompute
    dist_reset_flag <- FALSE;
    while (1)
    {
      bsamp <- sample (1:nrow (X), 2 * nrow (X), replace = TRUE, prob = kf_d$D / sum (kf_d$D)); # WARNING; twice the sample size
      while (1)
      {
        this_homer_cluster_method <- cluster_types[sample(length(cluster_types))[1]];
        this_homer_cluster_number <- cluster_numbers[sample(length(cluster_numbers))[1]];
        this_svm_kernel_type <- kernel_types[sample(length(kernel_types))[1]];
        if (feedback == TRUE) { cat (" (", substr (this_homer_cluster_method, 1, 3), ",", substr (this_svm_kernel_type, 1, 3), ") k = ", formatC (this_homer_cluster_number, 3)); }
        tmp_model <- try (homer (mldr_from_dataframe (cbind (X[bsamp,], Y[bsamp,]),  labelIndices=(ncol(X)+1):(ncol(X)+ncol(Y))), cores = 1, base.algorithm = method[2], method = this_homer_cluster_method, clusters = this_homer_cluster_number, scale = scale_flag, kernel = this_svm_kernel_type, cost = 10e+7, gamma = 10e-7), TRUE);
        
        if (class (tmp_model) != "try-error")
        {
          break;
        }
        this_homer_cluster_method <- "random";
        this_homer_cluster_number <- cluster_numbers[sample(length(cluster_numbers))[1]];
        cat ("Handling error\n");
        print (tmp_model);
      }
      this_pred        <- predict (tmp_model, X);
 
      if (save_in_disk == TRUE)
      {
        cat ("S");
        save (tmp_model, file = paste0 ("tmp_models/", "kfhe_ml", "_", paste0 (method, collapse = "_"), "_iter_", t, "_", unique_number, "_model.RData"));
        cat ("\bs ")
      }
      if (save_in_disk == FALSE)
      {
        kf_m$h_list[[t]] <- tmp_model;
      }
      rm (list = "tmp_model");
      gc ();
    
      
      # First transform the probabilities to classes, and then use these classes for the computation.
      if (blend_with_class == TRUE)
      {
        this_cls    <- factor (levels (Y)[apply (this_pred, 1, which.max)], levels = levels (Y));
        
        dummy_pred  <- matrix (0, nrow = nrow (this_pred), ncol = ncol (this_pred));
        colnames (dummy_pred) <- colnames (Y);
        for (i in 1:nrow(dummy_pred))
        {
          dummy_pred[i,this_cls[i]] <- 1;
        }
        
        this_pred <- dummy_pred;
      }
      # Now we have this_pred
      
      # For kf-m, calculate the measurement and the related error. This is a heuristic and can be computed in different ways.
      this_measurement <- (this_pred + kf_m$train_blend)/2;
      
      # Get the error for "this_measurement"
      unwt_comp_err  <- err_fun (this_measurement, Y, NULL); # Find the per datapoint error. This is a vector.
      uniwt_comp_err <- unwt_comp_err * (1/nrow (X));        # Get a uniformly weighted version. This is a weighted vector.
      

      # Also just the prediction error is computed for the learned model in 
      # this iteration, used to reset weights for datapoints
      this_pred_err       <- err_fun (this_pred, Y, NULL);
      wtd_this_pred_err   <- this_pred_err * kf_d$D;
      uniwt_this_pred_err <- this_pred_err * (1/nrow (X));
      
      this_m_err     <- uniwt_comp_err;
      this_d_err     <- uniwt_comp_err;
      
      
      # Measurement error for the model kf-m. This is a heuristic and can be computed in different ways.
      kf_m$V[t]          <- sum (this_m_err);# + log (length (levels (Y))-1);
      
      if (reset_dist == TRUE)
      {
        if ((sum (uniwt_this_pred_err) >= 0.5) && (dist_reset_flag == FALSE))
        {
          if (feedback == TRUE)
          {
            cat (", dreset = YES");
          }
          kf_d$D <- rep (1/nrow (X), nrow (X));
          kf_d$P_dist <- 1.0;
          
          dist_reset_flag <- TRUE;
          debug_info$reset_pt[t] <- dist_reset_flag;
          next;
        }
        else
        {
          if (dist_reset_flag == FALSE)
          {
            if (feedback == TRUE) 
            {
              cat (", dreset =  NO");
            }
          }
          break;
        }
      }
      else
      {
        if (feedback == TRUE) 
        {
          cat (", dreset =  NA");
        }
        break;
      }
    } # End of retry loop
    
    # Compute the Kalman gain for the kf-m
    kf_m$K[t]        <- kf_m$P[t] / (kf_m$P[t] + kf_m$V[t] + .Machine$double.xmin);
    # Blending the training predictions. This is not required for training, as we only need to store the kalman gains.
    # Update internal state for kf-m
    kf_m$train_blend <- kf_m$train_blend + kf_m$K[t] * (this_measurement - kf_m$train_blend);
    prev_blend_err   <- this_blend_err;
    
    # Update internal error for the kf-m
    P_t_pred <- kf_m$P[t] - kf_m$P[t] * kf_m$K[t];
    
    kf_m$P[t+1] <- P_t_pred;
    
    # Compute the actual error for the internal state
    debug_info$state_err[[t]]   <- err_fun(kf_m$train_blend, Y, NULL);
    this_blend_err              <- sum (err_fun(kf_m$train_blend, Y, NULL) * (1/nrow (X)));
    debug_info$blend_err_arr[t] <- this_blend_err;
    
    
    ##################################################################################################
    
    
    #############################################
    #### Weight Kalman filter (kf-w) section ####
    #############################################
    
    ############################
    #### Time update kf-w   ####
    ############################
#     # Based on present formulation, this is identity
    kf_d$D           <- kf_d$D;
    kf_d$P_dist      <- kf_d$P_dist + proc_noise;

    
    #################################
    #### Measurement update kf-w ####
    #################################
    # Measurement of state vector of the distribution kf-w
    dtemp           <- unwt_comp_err;
    dtemp <- dtemp + 1/nrow(X);
    kf_d$D_t_obs    <- kf_d$D * exp (dtemp);
    
    # Measurement error V_dist for the distribution kf-w
    kf_d$V_dist        <- sum (this_d_err);
 
    # Compute the Kalman gain for the distribution kf-w
    kf_d$K_dist[[t]]   <- kf_d$P_dist / (kf_d$P_dist + kf_d$V_dist + .Machine$double.xmin);

    # Update iternal state for the distribition kf-w
    kf_d$D             <- kf_d$D + kf_d$K_dist[[t]] * (kf_d$D_t_obs - kf_d$D); # Fishy
      
    # Update internal error for the distribution KF
    P_dist_old <- kf_d$P_dist;
    kf_d$P_dist        <- kf_d$P_dist - kf_d$P_dist * kf_d$K_dist[[t]];
    
    
    # Compute the change in the internal state actual error
    debug_info$blend_err_delta[t] <- prev_blend_err - this_blend_err;

    # Error feedback
    if (feedback == TRUE) 
    { 
      cat (", blend_err = ", formatC (this_blend_err, digits = print_prec, format = "f"), 
           "diff (", formatC (debug_info$blend_err_delta[t], digits = print_prec, format = "f"), ")", 
           "V_t = ", formatC (kf_m$V[t], digits = print_prec, format = "f"), ", ",
           "P_t = ", formatC (kf_m$P[t], digits = print_prec, format = "f"), ", ",
           "K_t = ", formatC (kf_m$K[t], digits = print_prec, format = "f"), ""); 
    }
  
    
    # Normalise distribution
    kf_d$D                   <- kf_d$D / sum (kf_d$D);
    if (feedback == TRUE)
    {
      cat (", ", "V_d = ", kf_d$V_dist, ", P_d = ", P_dist_old, "K_d = ", kf_d$K_dist[[t]]);
    }
    
    debug_info$D[[t]]        <- kf_d$D;
    
#     if ((early_stop == TRUE) && (this_blend_err < 10e-10))
    if ((early_stop == TRUE) && (kf_m$K[t] == 0)) # Similar in this case. Although practically it is better to stop on validation set.
    {
      max_iter <- t;
      break;
    }
    if (feedback == TRUE) { cat ("\n"); }
  }

  # Pack
  retmodel$kf_m             <- kf_m;
  retmodel$kf_d             <- kf_d;
  retmodel$debug_info       <- debug_info;
  retmodel$max_iter         <- max_iter;
  retmodel$cls_lvls         <- levels (Y);
  retmodel$blend_with_class <- blend_with_class;
  retmodel$method           <- method;
  retmodel$save_in_disk     <- save_in_disk;
  retmodel$unique_number    <- unique_number;
  
  # Stamp
  class (retmodel)  <- "kfhe_ml_m";
  
  # Send
  return (retmodel);
}



###########################################
####              Predict             #####
####################################################################################################
#### Function  : predict.kfhe_ml_m                                                              ####
#### Arguments :                                                                                #### 
####             model            : An object of type 'kfhe_m' returned by 'kfhe_train'         ####
####             X                : A matrix representing the datapoints to predict             ####
####             type             : If "prob", return the predicted per class scores            ####
####                                If "class", return the predicted class                      ####
####             feedback         : Print the states of the parameters while training           ####
####                                                                                            ####
#### Return    : Prediced scores or classes                                                     ####
####################################################################################################
predict.kfhe_ml_m <- function (model, X, type = "prob", feedback = TRUE)
{
  test_blend <- predict (model$kf_m$init_h, X);
  
  for (t in 1:model$max_iter)
  {
    if (model$save_in_disk == TRUE)
    {
      cat ("L");
      load (paste0 ("tmp_models/", "kfhe_ml", "_", paste0 (model$method, collapse = "_"), "_iter_", t, "_", model$unique_number, "_model.RData"));
      model$kf_m$h_list[[t]] <- tmp_model;
      cat ("\bl ");
    }
    
    if (feedback == TRUE) { cat ("\rPred iter = ", t, "      "); }
    this_pred <- predict (model$kf_m$h_list[[t]], X);
      
    # First transform the probabilities to classes, and then use these classes for the computation.
    if (model$blend_with_class == TRUE)
    {
      this_cls    <- factor (model$cls_lvls[apply (this_pred, 1, which.max)], levels = model$cls_lvls);
      
      dummy_pred  <- matrix (0, nrow = nrow (this_pred), ncol = ncol (this_pred));
      colnames (dummy_pred) <- model$cls_lvls;
      for (i in 1:nrow(dummy_pred))
      {
        dummy_pred[i,this_cls[i]] <- 1;
      }
      this_pred <- dummy_pred;
    }
    
    if (model$save_in_disk == TRUE)
    {
      model$kf_m$h_list[[t]] <- NULL;
      rm (list = "tmp_model");
      # TODO: Mechanism to remove the temporary models from disk, or just store them in /tmp
      gc ();
    }
    
    # This should be the same heuristic which is used during the training.
    this_measurement <- (test_blend + this_pred)/2;
    test_blend <- test_blend + model$kf_m$K[t] * (this_measurement - test_blend);
  }

  return (test_blend);
}

# TODO: Get rid of the unnecessary weight and othr things
err_fun <- function (pred, target, wt)
{
  if (!is.data.frame(target))
  {
    target         <- as.data.frame (target);
    one_hot_target <- model.matrix (~ ., data = target, contrasts.arg = lapply (target, contrasts, contrasts = FALSE));
    one_hot_target <- one_hot_target[,-1];
    target         <- one_hot_target;
  }
  else if (ncol (target) > 1)
  {
    
  }
  if (is.null (wt) == TRUE)
  {
    wt <- 1;
  }
  
  final_err <- wt * (rowMeans ((pred > 0.5) != target));
  
  return (final_err);
}

