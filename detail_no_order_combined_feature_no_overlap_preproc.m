function [feature_names,pos_feats,pos_labels,neg_feats,neg_labels] = detail_no_order_combined_feature_no_overlap_preproc(svm_preselect_method,kmer_length,kmer_distance,positiveFile,negativeFile,selected_kmers,revcom_kerms)
    %load positiveFile and negativeFile
	[feature_names,pos_feats,neg_feats] = getFeature(selected_kmers,revcom_kerms,positiveFile,negativeFile,svm_preselect_method,kmer_length,kmer_distance);
    
    input_size = size(pos_feats,2);
    disp(['Feature vector size: ',num2str(input_size)]);

	% TODO: change to measure feats matrixs
	posn = size(pos_feats,1);
	negn = size(neg_feats,1);

	%Generate labels
	pos_labels = ones(posn,1);
	neg_labels = ones(negn,1)*-1;
    %}
end

function  [kmercombi,duplication_mapping] = generateKmerCombiWithMapping(kmers,revcoms)
    kmercombi = {};
    duplication_mapping = [];
    index = 1;
    ss = 1;
    for c1 = 1:length(kmers)
        kmer1 = kmers(c1);
        revcom1 = revcoms(c1);
        for c2 = ss:length(kmers)
            kmer2 = kmers(c2);
            revcom2 = revcoms(c2);
            combination = generateAllCombination(kmer1,revcom1,kmer2,revcom2);
            
            kmercombi = horzcat(kmercombi,combination);
            duplication_mapping = horzcat(duplication_mapping,ones(1,length(combination))*index);
            
            index = index+1;
        end
        ss = ss+1;
    end
end

function combination = generateAllCombination(kmer1,rkmer1,kmer2,rkmer2)
    kmer1 = kmer1{1};
    rkmer1 = rkmer1{1};
    kmer2 = kmer2{1};
    rkmer2 = rkmer2{1};
    combination = {[kmer1,'+',kmer2],[kmer1,'+',rkmer2],[rkmer1,'+',kmer2],[rkmer1,'+',rkmer2]};
    combination = horzcat(combination,{[kmer2,'+',kmer1],[kmer2,'+',rkmer1],[rkmer2,'+',kmer1],[rkmer2,'+',rkmer1]});
    combination = unique(combination);
end

function [feats,selected_index] = extract_complex_feature(seqs,kmers_map,kmers_rev_map,kmer_length,combi_no,combis_map,combi_length,kmers_distance)
	input_size = length(seqs);

    no_extra_features = combi_no;

    %initialize feature space
    feats = zeros(input_size,no_extra_features);
    
    parfor i=1:input_size
        if mod(i,round(input_size/10)) == 0
            fprintf('%s%d/%d','Processing ',i,input_size);
        end
		
		%process the dna sequence
		seq = seqs{i};

        % initialize feature for current sample
        feat = zeros(1,no_extra_features);
        
        %get kmer index
        kmers_index_features = k_mers_index(seq,kmer_length,kmers_map);

        %cut the distance down if too long compare to dnalenght
        dna_length = length(kmers_index_features);
        if(kmers_distance > dna_length-kmer_length-1)
            c_kmers_distance = dna_length-kmer_length-1;
        else
            c_kmers_distance = kmers_distance;
        end

        for j = 1:(dna_length-1)
            if(kmers_index_features(j) ~= 0)
                %Get the kmer
                kmer1 = kmers_rev_map(kmers_index_features(j));

                %Decide second loop condition
                if(j+c_kmers_distance+kmer_length > dna_length)
                    c_length = dna_length;
                else
                    c_length = j+c_kmers_distance+kmer_length;
                end

                %Start second loop
                for k=j+kmer_length:c_length
                    if(kmers_index_features(k) ~= 0)
                        %Get the kmer
                        kmer2 = kmers_rev_map(kmers_index_features(k));
                        if(combi_length==2)
                            combi = [kmer1,'+',kmer2];
                            combiIndex = combis_map(combi);
                            feat(combiIndex) = feat(combiIndex)+1;
                        elseif(k < c_length && combi_length == 3)
                            % if select triplet combination
                            for l=k+1:c_length
                                if(kmers_index_features(l) ~= 0)
                                    %Get the kmer
                                    kmer3 = kmers_rev_map(kmers_index_features(l));
                                    combi = [kmer1,'+',kmer2,'+',kmer3];
                                    combiIndex = combis_map(combi);
                                    feat(combiIndex) = feat(combiIndex)+1;
                                end
                            end
                        end
                    end
                end
            end
        end

        %normalize feature vectors
        feat = feat/dna_length;
        %feat = feat/(dna_length*c_kmers_distance);
        %{
        if(sum(feat)>0)
            feat = feat/sum(feat);
        end
        %}

        feats(i,:) = feat;
    end
    
    %calculate non empty feature index
    count_total = sum(feats,1);
    selected_feature_bool = count_total > 0;
    feature_indexes = 1:length(combis_map);
    selected_index = feature_indexes(selected_feature_bool);
end

% returns all the k-mer in dna
function kmers_index_features = k_mers_index(dna,k,kmer_map)
	kmers_features = k_mers(dna,k);

	kmers_index_features = zeros(length(kmers_features)-k+1,1);
		
	for j=1:length(kmers_features)
		kmers_feature = kmers_features{j};
		if(kmer_map.isKey(kmers_feature))
			feat_index = kmer_map(kmers_feature);
			kmers_index_features(j) = feat_index;
		else
			kmers_index_features(j) = 0;
		end
	end
end

function [names,pos_feats,neg_feats] = getFeature(selected_kmers, revcom_kerms,pos_file,neg_file,svm_preselect_method,kmerlength,kmerdistance)
        kmers_maps = containers.Map(cat(2,selected_kmers,revcom_kerms), num2cell([1:length(selected_kmers),-(1:length(revcom_kerms))]));
        kmers_rev_map = containers.Map(num2cell([1:length(selected_kmers),-(1:length(revcom_kerms))]),cat(2,selected_kmers,revcom_kerms));

        %calculate combi length
        combi_length = str2double(svm_preselect_method(1));

        %Generate kmer combination
        disp('Generate kmer combination');
        [kmercombi,duplication_mapping] = generateKmerCombiWithMapping(selected_kmers,revcom_kerms);
        disp('Done!');

        %Select only unique combi
        %[kmercombi,unique_combi_index] = unique(kmercombi);
        %duplication_mapping = duplication_mapping(unique_combi_index);

        disp('Prepare feature mapping');
        %Select unique combi
        unique_combi_mapping = unique(duplication_mapping);
        unique_combi_no = length(unique_combi_mapping);

        combis_map = containers.Map(kmercombi, num2cell(duplication_mapping));
        disp('Done!');
        
		%preprocess the data
		disp('Extracting positive features...');
		[pos_feats,pos_feature_index] = preprocess(pos_file,kmers_maps,kmers_rev_map,kmerlength,unique_combi_no,combis_map,combi_length,kmerdistance);
        disp('Extracting negative features...');
        [neg_feats,neg_feature_index] = preprocess(neg_file,kmers_maps,kmers_rev_map,kmerlength,unique_combi_no,combis_map,combi_length,kmerdistance);
        disp('Done!');
        
        disp('Extract names');
        %calculate shared combi features
        for i=1:unique_combi_no
            for j=1:length(kmercombi)
                if(duplication_mapping(j) == i)
                    names{i} = kmercombi{j};
                    break;    
                end
            end
        end
        feature_index = union(pos_feature_index,neg_feature_index);
        pos_feats = pos_feats(:,feature_index);
        neg_feats = neg_feats(:,feature_index);
        names = names(feature_index);
        
        for i=1:length(names)
            names{i} = [names{i},9,revcom(names{i})];
        end
		disp('Done!');
end

function [feats,feat_indexes] = preprocess(file,kmers_maps,kmers_rev_map,kmerlength,combi_no,combis_map,combi_length,kmerdistance)
	file = ['sequence_files/',file,'.fa'];

	disp(['Loading input file ',file]);
	[seqs, ~] = loadDataFile(file);
	disp('Done!');

	%Extract features
	[feats,feat_indexes] = extract_complex_feature(seqs,kmers_maps,kmers_rev_map,kmerlength,combi_no,combis_map,combi_length,kmerdistance);
end

function [seqs,seqIds] = loadDataFile(dataFile)
	
	fid = fopen(dataFile);

	tline = fgetl(fid);

	i = 0;
	while ischar(tline)
		if tline(1) == '>'
            i = i + 1;
            
            seqIds{i} = tline(2:end);
            seqs{i} = [];
        else
            seqs{i} = [seqs{i},tline];
		end
		tline = fgetl(fid);
	end
	
	fclose(fid);
end

%get the revert combination of a kmer
function rc = revcom(kmer)
	rcbp = {{'A','T'}, {'G','C'}, {'C','G'}, {'T','A'},{'+','+'}};
	rc = [];
	for i = 1 : length(kmer)
		bp = kmer(i);
		for j=1:length(rcbp)
			if rcbp{j}{1} == bp
				rc = [rc,rcbp{j}{2}];
				break;
			end
		end
	end
	rc = rc(size(rc,2):-1:1);
end

function dic_map = createDicMap(dic)
	values = num2cell(1:length(dic));
	dic_map = containers.Map(dic, values);
end

%get the unique mapping table of all possiple kmers
function rcmapping = generateRevcompMappingTable(kmers)
	%create dictionary map of kmer
	idDicMap = createDicMap(kmers);

	l = 1;
	rcmapping = zeros(1,length(kmers));
	%for each kmer
	while l <= length(kmers)
		kmer = kmers{l};
		rc = revcom(kmer);

		%find the entry in dictionary
		kmerid = idDicMap(kmer);
        if(idDicMap.isKey(rc))
            rcid = idDicMap(rc);

            %check if the entry is found
            if rcid < kmerid
                rcmapping(l) = rcid;
            else
                rcmapping(l) = kmerid;
            end
        else
            rcmapping(l) = kmerid;
        end
		l = l+1;
	end
end

% returns all the k-mer in dna
function kmers_features = k_mers(dna, k)
	kmers_features = cell(1,length(dna)-k+1);
	for i = 1:(length(dna)-k+1)
		kmers_features{i} = upper(dna(i:i+k-1));
	end 
end