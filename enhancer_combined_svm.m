classdef enhancer_combined_svm < handle
	
	properties
		testing_predictions;
		testing_outputs;
		
		no_kmers;
		kmer_length;
		kmer_distance;
		preselect;
        svm_preselect_method;
		
		positive_file;
		negative_file;
        
        features_map;
        
        n_feature_names;
        n_pos_sample_inputs;
        n_pos_sample_outputs;
        n_neg_sample_inputs;
        n_neg_sample_outputs;
	end
	
	methods
		function a = enhancer_combined_svm(positive_file,negative_file,no_kmers,kmer_length,kmer_distance,preselect,svm_preselect_method,num_fold)
			a.no_kmers = no_kmers;
			a.kmer_length = kmer_length;
			a.kmer_distance = kmer_distance;
			a.preselect = preselect;
            a.svm_preselect_method = svm_preselect_method;
			
			a.positive_file = positive_file;
			a.negative_file = negative_file;
            
            a.features_map = 0;
			
            a.prepare_initial_data();
            
			a.n_fold(num_fold);
			
			%a.plot_roc()
		end
		
		function n_fold(a,n)
			pos_sample_size = size(a.n_pos_sample_inputs,1);
			neg_sample_size = size(a.n_neg_sample_inputs,1);

			pos_indices = crossvalind('Kfold',1:pos_sample_size,n);
			neg_indices = crossvalind('Kfold',1:neg_sample_size,n);
            
			%n=1;
			for i=1:n
				tic;
				disp(['Processing fold ',num2str(i),'/',num2str(n),' folds.']);
				
				pos_test = (pos_indices == i);
				neg_test = (neg_indices == i);

                %Extract features and split train and test data
                [f_name,testing_inputs,a.testing_outputs{i},training_inputs,training_outputs] = prepare_fold_data(a,pos_test,neg_test);				

				disp(['Size of training feartures vectors: ',num2str(size(training_inputs,1))]);
				disp(['Size of testing feartures vectors: ',num2str(size(testing_inputs,1))]);

				%Train the ensemble
				svm = a.train(training_inputs,training_outputs);
				
                clear training_inputs;
                clear training_outputs;
                
				toc;
				disp([num2str(i),'/',num2str(n),' folds processed.']);
                
                %Run test on this fold
                [a.testing_predictions{i},weight] = a.test(svm,testing_inputs);
                
                a.merge_feature_weight(f_name,weight);
            end
            
            clear a.pos_sample_inputs;
            clear a.pos_sample_outputs;
            clear a.neg_sample_inputs;
            clear a.neg_sample_outputs;
		end
		
		function [ROC_AUC,PR_AUC,feature_weight] = plot_roc(a,auc_only)
			n = length(a.testing_outputs);
            if(~auc_only)
                hold all;
                line_style = {'-','--','-.',':','-',':'};
                %color = {'k','k','k','k','k','r'};
            end
            
			ROC_AUC = 0;
            PR_AUC = 0;
            
			for i = 1:n
				% fit probabilities for scores
				%[X,Y,T,AUC] = perfcurve(ceil(a.testing_outputs{i}'/2),g_out{i}/2+0.5,1,'XVals',num2cell(0:0.001:1));
				[ROC_X,ROC_Y,~,R_AUC] = perfcurve(a.testing_outputs{i},a.testing_predictions{i},1);
				ROC_AUC = ROC_AUC+R_AUC;
                
                [PR_X,PR_Y,~,P_AUC] = perfcurve(a.testing_outputs{i},a.testing_predictions{i},1,'xCrit','reca','yCrit','prec');
                PR_AUC = PR_AUC+P_AUC;

                if(~auc_only)
                    % plot only if needed
                    p = plot(ROC_X,ROC_Y);
				
                    % Set DisplayNames for the lines for use by the legend
                    %disp(['Fold ',num2str(i),'/',num2str(n),' AUC:',num2str(AUC)]);
                    set(p, {'LineStyle'}, {line_style{i}});
                    %set(p, {'Color'},{color{i}});
                    set(p(1),'Displayname',['Fold ',num2str(i),'/',num2str(n),' ROC AUC:',num2str(R_AUC)]);
                    
                    % plot only if needed
                    p = plot(PR_X,PR_Y);
				
                    % Set DisplayNames for the lines for use by the legend
                    %disp(['Fold ',num2str(i),'/',num2str(n),' AUC:',num2str(AUC)]);
                    set(p, {'LineStyle'}, {line_style{i}});
                    %set(p, {'Color'},{color{i}});
                    set(p(1),'Displayname',['Fold ',num2str(i),'/',num2str(n),' PR AUC:',num2str(P_AUC)]);
                end
			end
			
			ROC_AUC = ROC_AUC/n;
            
            PR_AUC = PR_AUC/n;
            
            disp(['The average ROC AUC is ',num2str(ROC_AUC)]);
            
            disp(['The average PR AUC is ',num2str(PR_AUC)]);
            
            feature_weight = a.extract_feature_map(n);
			
            if(~auc_only)
                % Center a legend at the top of the graph
                legend('Location','SouthEast');
                xlabel('False positive rate and Recall'); ylabel('True positive rate and Precision');
                title('ROC and PR');

                hold off;
            end
        end
        
        function merge_feature_weight(a,f_names,weights)
            if(a.features_map ~= 0)
                %if exist then merge the new into the old
                for i=1:length(f_names)
                    f_name = f_names{i};
                    weight = weights(i);
                    if(a.features_map.isKey(f_name))
                        %if key exist
                        %update weight
                        a.features_map(f_name) = a.features_map(f_name)+weight;
                    else
                        %if key doesn't exist
                        %add key
                        a.features_map(f_name) = weight;
                    end
                end
            else 
                %if not exist initialize variable
                %create new map
                a.features_map = containers.Map(f_names,num2cell(weights));
            end
        end
        
        function [feature_weight] = extract_feature_map(a,n)
            feature_names = a.features_map.keys();
            
            weights = [];
            for i=1:length(feature_names)
                f_name = feature_names{i};
                weights(i) = a.features_map(f_name);
            end
            
            weights = weights/n;
            
            feature_weight = [feature_names',num2cell(weights')];

            feature_weight = sortrows(feature_weight,-2);
        end
		
		%Train SVM
        function  svm_struct = train(~,training_inputs,training_outputs)
			%a.svm_struct{i} = svmtrain(training_inputs,training_outputs,'boxconstraint',0.00001,'kernel_function','linear','options',statset('MaxIter',500000));
			%a.svm_struct{i} = svmtrain(training_inputs,training_outputs,'boxconstraint',0.0001,'kernel_function','linear','method','LS');
            svm_struct = svmtrain(training_inputs,training_outputs,'autoscale',false,'boxconstraint',100,'kernel_function','linear','method','SMO','options',statset('MaxIter',500000));
		end
		
		%Test SVM
		function [test_outputs,weight]=test(~,svm,test_inputs)
            %{
			shift = svm.ScaleData.shift;
			scale = svm.ScaleData.scaleFactor;
			test_inputs = bsxfun(@plus,test_inputs,shift);
			test_inputs = bsxfun(@times,test_inputs,scale);
            %}
			sv = svm.SupportVectors;
			alphaHat = svm.Alpha;
			bias = svm.Bias;
			kfun = svm.KernelFunction;
			kfunargs = svm.KernelFunctionArgs;
			test_outputs = kfun(sv,test_inputs,kfunargs{:})'*alphaHat(:) + bias;
			test_outputs = -test_outputs; % flip the sign to get the score for the +1 class
            
            weight  = -(alphaHat' * sv);
        end
        
        %Prepare single data, load normal data
        function prepare_initial_data(a)
            [a.n_feature_names,a.n_pos_sample_inputs,a.n_pos_sample_outputs,a.n_neg_sample_inputs,a.n_neg_sample_outputs] = data_preproc(a.kmer_length,a.positive_file,a.negative_file);
        end
		
		%Generate input data for each fold
		function [features_name,testing_inputs,testing_outputs,training_inputs,training_outputs] = prepare_fold_data(a,pos_test,neg_test)
            
            pos_train = ~pos_test;
            neg_train = ~neg_test;
            
            %preselect feature
            [selected_kmers,selected_revcoms] = preselect_feature(a,pos_train,neg_train);
            
            [c_feature_names,c_pos_sample_inputs,~,c_neg_sample_inputs,~] = detail_no_order_combined_feature_no_overlap_preproc(a.svm_preselect_method,a.kmer_length,a.kmer_distance,a.positive_file,a.negative_file,selected_kmers,selected_revcoms);
            
			pos_sample_inputs = [a.n_pos_sample_inputs,c_pos_sample_inputs];
			pos_sample_outputs = a.n_pos_sample_outputs;
			neg_sample_inputs = [a.n_neg_sample_inputs,c_neg_sample_inputs];
			neg_sample_outputs = a.n_neg_sample_outputs;
            features_name = [a.n_feature_names,c_feature_names];
			%{}
			
			
			%[a.pos_sample_inputs,a.pos_sample_outputs,a.neg_sample_inputs,a.neg_sample_outputs] = combined_feature_preproc(a.preselect,a.no_kmers,a.kmer_length,a.kmer_distance,a.positive_file,a.negative_file);

			%[a.pos_sample_inputs,a.pos_sample_outputs,a.neg_sample_inputs,a.neg_sample_outputs] = data_preproc(a.kmer_length,a.positive_file,a.negative_file);
			
			%normalize
			%{}
			for i=1:size(pos_sample_inputs,1)
				pos_sample_inputs(i,:) = pos_sample_inputs(i,:)/sum(pos_sample_inputs(i,:));
			end
			
			for i=1:size(neg_sample_inputs,1)
				neg_sample_inputs(i,:) = neg_sample_inputs(i,:)/sum(neg_sample_inputs(i,:));
			end
			%{}
            
            
            %split train and test data
            pos_train = ~pos_test;
            neg_train = ~neg_test;
            
            pos_test_set_inputs = pos_sample_inputs(pos_test,:);
            neg_test_set_inputs = neg_sample_inputs(neg_test,:);

            testing_inputs = [pos_test_set_inputs;neg_test_set_inputs];
            clear pos_test_set_inputs;
            clear neg_test_set_inputs;

            pos_train_set_inputs = pos_sample_inputs(pos_train,:);
            neg_train_set_inputs = neg_sample_inputs(neg_train,:);

            training_inputs = [pos_train_set_inputs;neg_train_set_inputs];
            clear pos_train_set_inputs;
            clear neg_train_set_inputs;

            pos_test_set_outputs = pos_sample_outputs(pos_test,:);
            neg_test_set_outputs = neg_sample_outputs(neg_test,:);
            pos_train_set_outputs = pos_sample_outputs(pos_train,:);
            neg_train_set_outputs = neg_sample_outputs(neg_train,:);

            testing_outputs = [pos_test_set_outputs;neg_test_set_outputs];

            training_outputs = [pos_train_set_outputs;neg_train_set_outputs];
        end
        
        
        %depend on selection method, preselect k-mers
        function [selected_kmers,revcom_kmers] = preselect_feature(a,pos_train,neg_train)
            if(strcmp(a.preselect, 'adaboost'))
                %Preselect using adaboost
                disp('Pre select features...');
                
                adaboost_struct = pre_select_enhancer_adaboost(a.positive_file,a.negative_file,a.kmer_length,50,2,a.no_kmers,pos_train,neg_train);
                selected_kmers = adaboost_struct.get_selected_kmers();
                
                %split selected kmers and revcom
                [selected_kmers,revcom_kmers] = a.split(selected_kmers);
                
                disp('Done!');
            elseif(strcmp(a.preselect, 'svm'))
                %Preselect using svm
                disp('Pre select features...');
                
                svm_struct = pre_select_enhancer_svm(a.positive_file,a.negative_file,a.kmer_length,pos_train,neg_train);
                [~,~,feature_weight] = svm_struct.plot_roc(true);
                selected_kmers = feature_weight(:,1);
                
                %split selected kmers and revcom
                [selected_kmers,revcom_kmers] = a.split(selected_kmers);

                disp('Done!');
            else
                %Preselect using other preselection methods
                disp('Pre select features...');
                
                struct = pre_select_enhancer(a.positive_file,a.negative_file,a.kmer_length,a.preselect,a.no_kmers,pos_train,neg_train);
                selected_kmers = struct.selected_features;
                
                %split selected kmers and revcom
                [selected_kmers,revcom_kmers] = a.split(selected_kmers);
                
                disp('Done!');
            end
            
            %Select top weight if +
            %Select both top and bottom weight if +-
            %Select bottom weight if -
            if(~isempty(strfind(a.svm_preselect_method,'+-')))
                if(mod(a.no_kmers,2)==0)
                    half = a.no_kmers/2;
                    l = length(selected_kmers);
                    selected_kmers = selected_kmers([1:half,(l-half+1):l]);
                    revcom_kmers = revcom_kmers([1:half,(l-half+1):l]);
                else
                    error(['Do not support odd selected_no with pre select method ',a.svm_preselect_method]);
                end
            elseif(strcmp(a.svm_preselect_method(end-1),'+'))
                selected_kmers = selected_kmers(1:a.no_kmers);
                revcom_kmers = revcom_kmers(1:a.no_kmers);
            elseif(strcmp(a.svm_preselect_method(end-1),'-'))
                l = length(selected_kmers);
                selected_kmers = selected_kmers((l-a.no_kmers+1):l);
                revcom_kmers = revcom_kmers((l-a.no_kmers+1):l);
            end
        end
        
        %depend on selection method, preselect k-mers
        function [s_kmers,s_revcoms] = split(a,selected_kmers) 
            for i=1:length(selected_kmers)
                tline = selected_kmers{i};
                [kmer,reman] = strtok(tline,9);
                [recom,~] = strtok(reman(2:end),9);

                s_kmers{i} = kmer;
                s_revcoms{i} =  recom;
            end
        end
		
	end
	
end

