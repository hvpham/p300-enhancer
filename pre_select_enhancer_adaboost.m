classdef pre_select_enhancer_adaboost < handle
	
	properties
		ens;

        kmers;
        
		pos_sample_inputs;
		pos_sample_outputs;
		neg_sample_inputs;
		neg_sample_outputs;
        
        testing_inputs;
        testing_outputs;
        
		kmer_length;
        initial_size;
        resume_step;
        
        no_selected;
        
		positive_file;
		negative_file;
        
        pos_use;
        neg_use;
	end
	
	methods
		
		function a = pre_select_enhancer_adaboost(positive_file,negative_file,kmer_length,initial_ensemble_size,resume_step,no_selected,pos_use,neg_use)
			a.kmer_length = kmer_length;
            
            a.kmer_length = kmer_length;
            a.initial_size = initial_ensemble_size;
            a.resume_step = resume_step;
            
            a.no_selected = no_selected;
            
			a.positive_file = positive_file;
			a.negative_file = negative_file;
            
            a.pos_use = pos_use;
            a.neg_use = neg_use;
			
			%Load data
			a.get_data();
			
            %Split data for testting
            [train_inputs,train_outputs] = a.split_data();
            
            %Train the ensemble
            a.train(train_inputs,train_outputs);
            
            a.resume(a.resume_step)
        end
        
        function resume(a,step)
            no_seleted_kmers = a.count_selected_kmers();
            
            resume = true;
            i=1;
            
            while(resume)
                if(no_seleted_kmers >= a.no_selected)
                    disp('Done optimizing');
                    break;
                end
                                
                a.ens = a.ens.resume(step);
                
                no_seleted_kmers = a.count_selected_kmers();
                
                disp(['Finish round ',num2str(i)]);

                i=i+1;
            end
        end
        
        %Split data for testting
        function [train_inputs,train_outputs] = split_data(a)
            pos_sample_size = size(a.pos_sample_inputs,1);
			neg_sample_size = size(a.neg_sample_inputs,1);

			pos_indices = crossvalind('Kfold',1:pos_sample_size,5);
			neg_indices = crossvalind('Kfold',1:neg_sample_size,5);

            pos_test = (pos_indices == 1); pos_train = ~pos_test;
            neg_test = (neg_indices == 1); neg_train = ~neg_test;

            pos_test_set_inputs = a.pos_sample_inputs(pos_test,:);
            neg_test_set_inputs = a.neg_sample_inputs(neg_test,:);
            pos_train_set_inputs = a.pos_sample_inputs(pos_train,:);
            neg_train_set_inputs = a.neg_sample_inputs(neg_train,:);

            pos_test_set_outputs = a.pos_sample_outputs(pos_test,:);
            neg_test_set_outputs = a.neg_sample_outputs(neg_test,:);
            pos_train_set_outputs = a.pos_sample_outputs(pos_train,:);
            neg_train_set_outputs = a.neg_sample_outputs(neg_train,:);

            a.testing_inputs = [pos_test_set_inputs;neg_test_set_inputs];
            a.testing_outputs = [pos_test_set_outputs;neg_test_set_outputs];
            train_inputs = [pos_train_set_inputs;neg_train_set_inputs];
            train_outputs = [pos_train_set_outputs;neg_train_set_outputs];
            
        end
        
        %Check if there are enought selected kmer
		function no_selected_kmers = count_selected_kmers(a)
            trained_tree = a.ens.Trained;
            kmers_index = zeros(length(trained_tree),1);
            
            for i=1:length(trained_tree)
                cut_var = trained_tree{i}.CutVar;
                kmers_index(i) = str2num(cut_var{1}(2:end));
            end
            
            no_selected_kmers = length(unique(kmers_index(1:i)));
            
			disp([num2str(no_selected_kmers),' selected kmers']);
		end
	
		%Train SVM
		function train(a,training_inputs,training_outputs)
            a.ens = fitensemble(training_inputs,training_outputs,'AdaBoostM1',a.initial_size,'Tree');
        end
        
        %Test Adaboost
        function outputs = test(a,testing_inputs)
            [~,outputs] = a.ens.predict(testing_inputs);
        end
        
        %Get selected kmers
        function [selected_kmers] = get_selected_kmers(a)
            trained_tree = a.ens.Trained;
            kmers_index = zeros(length(trained_tree),1);
            for i=1:length(trained_tree)
                cut_var = trained_tree{i}.CutVar;
                kmers_index(i) = str2num(cut_var{1}(2:end));
                
                if(length(unique(kmers_index(1:i))) == a.no_selected)
                    kmers_index = unique(kmers_index(1:i));
                    break;
                end
            end
            
            selected_kmers = a.kmers(kmers_index);
            
            %{
            w =  abs(w);

            [w,wI] = sort(w);

            selected_kmers = reIndex(selected_kmers,wI);
            selected_kmers = reIndex(selected_kmers,[1:selected_no]);

            revcom_kmers = reIndex(revcom_kmers,wI);
            revcom_kmers = reIndex(revcom_kmers,[1:selected_no]);
            %}

            %{
            kmers_length = length(a.kmers);

            selected_kmers = a.kmers([1:(selected_no/2),(kmers_length-selected_no/2+1):kmers_length]);

            revcom_kmers = a.r_kmers([1:selected_no/2,kmers_length-selected_no/2+1:kmers_length]);
            %}

            %{
            selected_kmers = reIndex(selected_kmers,[1:(selected_no)]);

            revcom_kmers = reIndex(revcom_kmers,[1:selected_no]);
            %}
        end
        
		%Generate some data from the required function.
		function get_data(a)
            [a.kmers,a.pos_sample_inputs,a.pos_sample_outputs,a.neg_sample_inputs,a.neg_sample_outputs] = data_preproc(a.kmer_length,a.positive_file,a.negative_file);
            
            %use only the portion allowed
            a.pos_sample_inputs = a.pos_sample_inputs(a.pos_use,:);
            a.pos_sample_outputs = a.pos_sample_outputs(a.pos_use,:);
            a.neg_sample_inputs = a.neg_sample_inputs(a.neg_use,:);
            a.neg_sample_outputs = a.neg_sample_outputs(a.neg_use,:);
            
            %normalize
            for i=1:size(a.pos_sample_inputs,1)
                a.pos_sample_inputs(i,:) = a.pos_sample_inputs(i,:)/sum(a.pos_sample_inputs(i,:));
            end
            
            for i=1:size(a.neg_sample_inputs,1)
                a.neg_sample_inputs(i,:) = a.neg_sample_inputs(i,:)/sum(a.neg_sample_inputs(i,:));
            end
        end
	end
	
end

