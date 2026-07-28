knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# install.packages("pak")
# pak::pak("bennett-stolze/NestNetR")

library(NestNetR)

ID <- "BB959"
Species <- "RuddyTurnstone" 
wd <- "data"
dir.raw <- file.path(wd, "RawData", Species)

raw_light <- read_light(file.path(wd, "RawData", Species, paste0(ID, ".lux")))
head(raw_light)
str(raw_light)

raw_deg <- read_deg(file.path(wd, "RawData", Species, paste0(ID, ".deg")))
head(raw_deg)
str(raw_deg)

tm.breeding <- set_breeding_period(dir.raw, raw_light, raw_deg, ID, auto = TRUE)
tm.breeding

# tm.breeding <- set_breeding_period(dir.raw, raw_light, raw_deg, ID, auto = FALSE, gr.Device = "x11")
# tm.breeding

segment_days <- 1
breeding_data <- preprocessing(ID, raw_light, raw_deg, tm.breeding, segment_days = segment_days)

# library(keras3)
# library(caret)
#
# # Define window length (number of days)
# segment_days <- 1
#
# # Define class weights
# weight_brooding <- 2.5
# weight_incubation <- 2
# weight_random <- 0.8

# # gather data from all available geolocator records of your species
# breeding_data_list <- create_breeding_data_list(dir.raw, segment_days = segment_days)
#
# # Create training data
# preclassified_data <- create_trainingdata(breeding_data_list,
#                                           segment_days = segment_days,
#                                           dir.raw = dir.raw,
#                                           gr.Device = "x11")
# 
# # Partition preclassified data into training- (80%) and test-data (20%)
# classes <- sapply(preclassified_data, `[[`, "Class") # extract classes
# partition <- caret::createDataPartition(classes, p = 0.8, list = FALSE, times = 1) # 80% of each class goes to training, 20% to test data
# training_data <- preclassified_data[partition]
# test_data <- preclassified_data[-partition]

# vars <- c("Light", "Tmin", "Tmax")
# window_length <- segment_days * 288
#
# make_array <- function(data, vars, window_length) {
#   channels <- lapply(vars, function(var) {
#     do.call(rbind, lapply(data, function(x) as.numeric(x[[var]])[seq_len(window_length)]))
#   })
#   do.call(abind::abind, c(channels, along = 3))
# }
# 
# x_breeding <- make_array(breeding_data, vars, window_length)
# x_train <- make_array(training_data, vars, window_length)
# x_test <- make_array(test_data, vars, window_length)

# model_name <- "custom"
# # link directory
# dir.model <- file.path(wd, "Model", model_name)
# 
# # Create directory inside your data folder for storing the custom model
# dir.create(dir.model, recursive = TRUE, showWarnings = FALSE)
# 
# # Create k folds
# class_levels <- c("brooding", "incubation", "random")
# train_classes <- factor(sapply(training_data, `[[`, "Class"), levels = class_levels)
# folds <- caret::createFolds(train_classes, k = 5, returnTrain = TRUE)
# 
# # K-fold cross-validation
# cv_results <- lapply(seq_along(folds), function(i) {
#   cat("Processing Fold", i, "of", length(folds), "\n")
# 
#   # Get indices for this fold
#   train_idx <- folds[[i]]
#   val_idx <- setdiff(seq_len(length(training_data)), train_idx)
# 
#   # Split the data
#   x_train_fold <- x_train[train_idx, , , drop = FALSE]
#   x_val_fold <- x_train[val_idx, , , drop = FALSE]
# 
#   # Prepare labels (one-hot encode)
#   y_train_fold <- to_categorical(as.integer(train_classes[train_idx]) - 1,
#                                  num_classes = length(class_levels))
#   y_val_fold <- to_categorical(as.integer(train_classes[val_idx]) - 1,
#                                num_classes = length(class_levels))
# 
# 
#   # Define the model
#   inputs <- layer_input(shape = c(window_length, dim(x_train)[3]))
#   outputs <- inputs %>%
#     layer_conv_1d(filters = 400, kernel_size = 12, activation = 'relu') %>%
#     layer_average_pooling_1d(pool_size = 12, padding = "same") %>%
#     layer_conv_1d(filters = 288, kernel_size = 16, activation = 'relu') %>%
#     layer_average_pooling_1d(pool_size = 6, padding = "same") %>%
#     layer_flatten() %>%
#     layer_dropout(rate = 0.2) %>%
#     layer_dense(units = 160, activation = 'relu') %>% #190 or 128
#     layer_dense(units = ncol(y_train_fold), activation = 'softmax')
#   model <- keras_model(inputs = inputs,
#                        outputs = outputs)
# 
#   # Compile the model
#   model %>% compile(
#     loss = 'categorical_crossentropy',
#     optimizer = optimizer_adam(
#       learning_rate = 0.001,
#       beta_1 = 0.9,
#       beta_2 = 0.999,
#       epsilon = 1e-7
#     ),
#     metrics = c('categorical_accuracy', 'Precision', 'Recall')
#   )
# 
#   # Train model
#   history <- model %>% keras3::fit(
#     x_train_fold, y_train_fold,
#     validation_data = list(x_val_fold, y_val_fold),
#     epochs = 50,
#     batch_size = 32,
#     verbose = 0,
#     class_weight = list(`0` = weight_brooding,    # brooding
#                         `1` = weight_incubation,  # incubation
#                         `2` = weight_random),     # random
#     callbacks = list(
#       callback_early_stopping(monitor = "val_loss", patience = 8, restore_best_weights = TRUE),
#       callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 4, min_lr = 1e-6),
#       callback_model_checkpoint(filepath = file.path(dir.model, paste0("cv_fold_", i, ".keras")), save_best_only = TRUE, monitor = "val_loss")
#       )
#   )
# 
#   # Evaluate on validation fold
#   val_metrics <- model %>% evaluate(x_val_fold, y_val_fold, verbose = 0)
# 
#   val_accuracy  <- as.numeric(val_metrics["categorical_accuracy"])
#   val_precision <- as.numeric(val_metrics["Precision"])
#   val_recall    <- as.numeric(val_metrics["Recall"])
#   val_f1 <- 2 * (val_precision * val_recall) / (val_precision + val_recall)
#
#   cat("Fold", i, "- Val F1:", round(val_f1, 2), "\n")
#   keras3::clear_session() # clear between folds
#
#   data.frame(
#     fold = i,
#     val_accuracy,
#     val_precision,
#     val_recall,
#     val_f1
#   )
# }) |>
#   dplyr::bind_rows()
# 
# # Summary of cross-validation results
# cat("\n=== Cross-Validation Results ===\n")
# print(cv_results)
# cat("Mean Validation F1:", round(mean(cv_results$val_f1, na.rm = TRUE), 4), "±", round(sd(cv_results$val_f1, na.rm = TRUE), 4), "\n")
# # Precision & Recall = averages over all classes
# 
# # Examine CV results for stability and performance
# cat("CV Performance Summary:\n")
# print(summary(cv_results))
# 
# # Check for high variance across folds (indicates instability)
# if (sd(cv_results$val_accuracy) > 0.1) {
#   cat("Warning: High variance across folds. Consider:\n")
#   cat("- More data\n- Different architecture\n- Better regularization\n")
# }

# # indicate file path for final model
# path_final_model <- file.path(dir.model, "final_model.keras")
# 
# # Convert labels to numeric and one-hot encode
# classes_training <- sapply(training_data, `[[`, "Class")
# y_train <- to_categorical(as.integer(factor(classes_training, levels = class_levels)) - 1,
#                           num_classes = length(class_levels))
# 
# # Load the best performing model from Cross-Validation
# model <- keras3::load_model(file.path(dir.model, paste0("cv_fold_", cv_results$fold[which.max(cv_results$val_f1)], ".keras")))
# summary(model)
# # Compile the model
# model %>% compile(
#   loss = 'categorical_crossentropy',
#   optimizer = optimizer_adam(
#     learning_rate = 0.001,
#     beta_1 = 0.9,
#     beta_2 = 0.999,
#     epsilon = 1e-7),
#   metrics = c('categorical_accuracy', 'Precision', 'Recall')
#   )
# 
# history <- model %>% keras3::fit(
#   x_train, y_train,
#   epochs = 50,
#   batch_size = 64,
#   class_weight = list(`0` = weight_brooding,   # brooding
#                     `1` = weight_incubation,   # incubation
#                     `2` = weight_random),      # random
#   callbacks = list(
#       callback_early_stopping(monitor = "loss", patience = 3, restore_best_weights = TRUE),
#       callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.5, patience = 4, min_lr = 1e-4),
#       callback_model_checkpoint(filepath = path_final_model, save_best_only = TRUE, monitor = "loss")
#   )
# )

# # Make predictions on test set
# classes_test <- sapply(test_data, `[[`, "Class")
# y_test <- to_categorical(as.integer(factor(classes_test, levels = class_levels)) - 1,
#                          num_classes = length(class_levels))
# 
# test_predictions <- model %>% predict(x_test)
# test_metrics <- model %>% evaluate(x_test, y_test)
# 
# # Detailed evaluation
# pred_classes <- apply(test_predictions, 1, which.max) - 1
# true_classes <- apply(y_test, 1, which.max) - 1
# 
# # Confusion matrix
# cm <- confusionMatrix(
#   factor(pred_classes, levels = 0:2, labels = c("brooding", "incubation", "random")),
#   factor(true_classes, levels = 0:2, labels = c("brooding", "incubation", "random"))
# )

classified_breeding <- classify_breeding_behaviour(breeding_data)
classified_breeding

# classified_breeding <- classify_breeding_behaviour(breeding_data, model = path_final_model)
# classified_breeding

# # Define species and directories
# Species <- "RuddyTurnstone"
# wd <- "data"
# dir.raw <- file.path(wd, "RawData", Species)
# 
# # gather data from all available geolocator records of your species
# breeding_data_list <- create_breeding_data_list(dir.raw, segment_days = 1)
# 
# # remove empty lists
# breeding_data_list <- breeding_data_list[lengths(breeding_data_list) > 0]
# 
# classified_breeding_df <- lapply(breeding_data_list, classify_breeding_behaviour) |>
#   dplyr::bind_rows()
# 
# # Save results to CSV
# utils::write.csv(
#   classified_breeding_df,
#   file.path(wd, paste0("classified_breeding_results_", Species, ".csv")),
#   row.names = FALSE
# )
