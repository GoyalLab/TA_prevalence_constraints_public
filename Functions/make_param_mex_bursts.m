% original version written by Lea Schuh and Yogesh Goyal, 2020
% Schuh et al., Cell Systems, 2020

% Create a C file running the Gillespie simulations as well as a
% parameter.m file summarizing the parameters and propensities. 

function uniquetok = make_param_mex_bursts(inputfile, outputfile, n_species, iruns)

% First thing is to read in inputfile
fid = fopen(inputfile);

% finc = fopen('/storage/scratch/users/lea.schuh/RajNetworkModeling/basicgillespie.c');
finc = fopen('basicgillespie.c');

% The first line is the number of Gillespie steps to do
gillstepsline = fgetl(fid);
maxgillespiesteps = str2double(gillstepsline);

clear inlines lines rate ratenames ratelines

inlines = cell(1,4*n_species);
onlines = cell(1,1);
b_lines = cell(1,n_species);
k_nlines = cell(1, n_species);
nlines = cell(1,n_species);
initiallines = cell(1,n_species);
parmat_line = cell(1,n_species*n_species);
parmat_linnum = zeros(1,n_species*n_species);

i=1;
while 1 %Rates
    tline = fgetl(fid); 
    while isempty(tline)
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Production difference'), break, end
    inlines{i}=tline;
    i=i+1;
end

i=1;
while 1 %Difference in production 
    tline = fgetl(fid);
    while isempty(tline)
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Basal values'), break, end
    onlines{i}=tline;
    i=i+1;
end

i=1; %Basal rates
while 1
    tline = fgetl(fid);
    while isempty(tline) 
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Hill function k'), break, end
    b_lines{i}=tline;
    i=i+1;
end

i=1;
while 1  %Hill function k
    tline = fgetl(fid);
    while isempty(tline)
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Hill function n'), break, end
    k_nlines{i} = tline;
    i=i+1;
end

i=1;
while 1  %Hill function n
    tline = fgetl(fid);
    while isempty(tline)
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Initial values'), break, end
    nlines{i} = tline;
    i=i+1;
end

i=1;
while 1  %Initial values
    tline = fgetl(fid);
    while isempty(tline) 
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    if strcmp(tline,'Network'), break, end
    initiallines{i} = tline;
    i=i+1;
end

i=1;
while 1  %Network structure
    tline = fgetl(fid);
    while isempty(tline) 
        tline = fgetl(fid);
    end
    if ~ischar(tline), break, end
    parmat_line{i} = tline;
    i=i+1;
end

%Retrieve the network matrix
for ipar_num = 1:length(parmat_line)
    parmat_linnum(ipar_num) = str2double(parmat_line{ipar_num});
end
parmat = zeros(n_species);
for ipar = 1:n_species
    parmat(ipar,:) = parmat_linnum(n_species*ipar-(n_species-1):n_species*ipar);
end

%Extract the rates from the beginning of the lines
ratelines = cell(1,4*n_species);            %Rates
lines = cell(1,4*n_species);
ratenames = cell(1,4*n_species);
rate= zeros(1,4*n_species);

on_ratelines = cell(1,n_species);           %Difference in production arte due to on-state of burst
on_ratenames = cell(1,n_species);
on_rate = zeros(1,n_species);

b_ratelines = cell(1,n_species);            %Basal values
b_ratenames = cell(1,n_species);
b_rate = zeros(1,n_species);

k_ratelines = cell(1,n_species);            %Dissociation constant k (Hill function)
k_ratenames = cell(1,n_species);
k_rate = zeros(1,n_species);

n_ratelines = cell(1,n_species);            %Hill coefficient k (Hill function)
n_ratenames = cell(1,n_species);
n_rate = zeros(1,n_species);

j=1;

for i = 1:length(inlines)
    clear R K
    [ratelin,lin] = strread(inlines{i},'%s%s','delimiter',':');
    ratelines{i} = strrep(ratelin{1},',',' ');  %Replace comma with space
    lines{j} = lin{1};
    [R,K] = strread(ratelines{i},'%s%f','delimiter','=');
    ratenames{j} = R{1};
    rate(j) = K(1);
    if length(R)==2  %If it's bidirectional...
        % First, split the line into the lhs and the rhs:
        j = j+1;  %Now we add the rates and the reverse rxn:
        ratenames{j} = R{2};
        rate(j) = K(2);
        [LHS,RHS] = strread(lin{1},'%s%s','delimiter','=');
        lines{j} = [RHS{1},'=',LHS{1}];
    end
    j = j+1;
end

j=1;
for i = 1:length(onlines) %Rates
    clear R K
    on_ratelines{i} = strrep(onlines{i},',',' '); 
    [R,K] = strread(on_ratelines{i},'%s%f','delimiter','=');
    on_ratenames{j} = R{1};
    on_rate(j) = K(1);
    j = j+1;
end

j=1;
for i = 1:length(b_lines) %Basal rates
    clear R K
    b_ratelines{i} = strrep(b_lines{i},',',' ');
    [R,K] = strread(b_ratelines{i},'%s%f','delimiter','=');
    b_ratenames{j} = R{1};
    b_rate(j) = K(1);
    j = j+1;
end

j=1;
for i = 1:length(k_nlines) %Hill function k
    clear R K
    k_ratelines{i} = strrep(k_nlines{i},',',' '); 
    [R,K] = strread(k_ratelines{i},'%s%f','delimiter','=');
    k_ratenames{j} = R{1};
    k_rate(j) = K(1);
    j = j+1;
end

j=1;
for i = 1:length(nlines) %Hill function n
    clear R K
    n_ratelines{i} = strrep(nlines{i},',',' ');
    [R,K] = strread(n_ratelines{i},'%s%f','delimiter','=');
    n_ratenames{j} = R{1};
    n_rate(j) = K(1);
    j = j+1;
end

fclose(fid);

numrxns = length(lines);
currrxn = 1;

clear alllhs allrhs alltok rxntok

for currrxn = 1:numrxns
    
    clear lhstok rhstok;
    
    % First, split the line into the lhs and the rhs:
    [LHS,RHS] = strread(lines{currrxn},'%s%s','delimiter','=');
    
    % Now get all tokens from the LHS...
    if ~isempty(LHS)
        lhstok = strtrim(strread(LHS{1},'%s','delimiter','+'));
        rxntok{currrxn}.lhs = lhstok; %remove whitespace
    else
        lhstok = {};
        rxntok{currrxn}.lhs = {};
    end
    % Now get all tokens from the LHS...
    if ~isempty(RHS)
        rhstok = strtrim(strread(RHS{1},'%s','delimiter','+'));
        rxntok{currrxn}.rhs = rhstok;
        %remove whitespace
    else
        rhstok = {};
        rxntok{currrxn}.rhs = {};
    end
    alltok{currrxn} = [lhstok',rhstok'];
end

% Make a list of all the species we have
tok = [alltok{:}];
uniquetok = unique(tok,'stable');
nspecies = length(uniquetok);


if iruns == 1
    
%     foutc = fopen(['/storage/scratch/users/lea.schuh/RajNetworkModeling/' outputfile 'histomex.c'],'w');
    foutc = fopen([outputfile 'histomex.c'],'w');
    
    %Create C file 
    lin = 'junk';
    while ~strcmp(lin,'//INCLUDE NSPECIES HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end
    fprintf(foutc,'#define NSPECIES %d\n',nspecies);

    while ~strcmp(lin,'//INCLUDE NUMRXNS HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end
    fprintf(foutc,'#define NUMRXNS %d\n',numrxns);

    while ~strcmp(lin,'//INSERT ALL VARIABLE DECLARATIONS HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end

    fprintf(foutc,'  long %s',uniquetok{1});
    for i = 2:nspecies
        fprintf(foutc,',%s',uniquetok{i});
    end
    fprintf(foutc,';\n');

    fprintf(foutc,'  double %s',ratenames{1});

    for i = 2:length(ratenames)
        fprintf(foutc,',%s',ratenames{i});
    end

    for i = 1:length(on_ratenames)
        fprintf(foutc,',%s',on_ratenames{i});
    end

    for i = 1:length(b_ratenames)
        fprintf(foutc,',%s',b_ratenames{i});
    end

    for i = 1:length(k_ratenames)
        fprintf(foutc,',%s',k_ratenames{i});
    end

    for i = 1:length(n_ratenames)
        fprintf(foutc,',%s',n_ratenames{i});
    end

    fprintf(foutc,';\n');

    while ~strcmp(lin,'//UNPACK ALL SPECIES HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end

    for i = 1:nspecies
        fprintf(foutc,'  %s = (long)species[%d];\n',uniquetok{i},i-1);
    end

    while ~strcmp(lin,'//UNPACK ALL RATES HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end

    for i = 1:length(ratenames)
        fprintf(foutc,'  %s = rates[%d];\n',ratenames{i},i-1);
    end
    for d = 1:length(on_ratenames)
        fprintf(foutc,'  %s = rates[%d];\n',on_ratenames{d},i+d-1);
    end
    for l = 1:length(b_ratenames)
        fprintf(foutc,'  %s = rates[%d];\n',b_ratenames{l},i+d+l-1);
    end
    for j = 1:length(k_ratenames)
        fprintf(foutc,'  %s = rates[%d];\n',k_ratenames{j},i+d+l+j-1);
    end
    for k = 1:length(n_ratenames)
        fprintf(foutc,'  %s = rates[%d];\n',n_ratenames{k},i+d+l+j+k-1);
    end

    clear lhsrxnlist

    %Define the stoichiometry matrix.
    lhsstoichiometry = zeros(nspecies,numrxns);
    rhsstoichiometry = zeros(nspecies,numrxns);
    for i = 1:nspecies
        for j = 1:numrxns
            lhsstoichiometry(i,j) = sum( ismember(rxntok{j}.lhs, uniquetok(i)) );
            rhsstoichiometry(i,j) = sum( ismember(rxntok{j}.rhs, uniquetok(i)) );
        end
    end

    %*****************************
    %Timestep sampling
    %Save data
    while ~strcmp(lin,'//INSERT SAVE HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end
    %Now let's generate code for saving the data:
    for j=1:nspecies
        fprintf(foutc,'species_out[savecount*NSPECIES+%d] = %s;\n',j-1,uniquetok{j});
    end
    while ~strcmp(lin,'//INSERT IF STATEMENT HERE')
        lin = fgetl(finc);
        fprintf(foutc,'%s\n',lin);
    end

    for i = 1:numrxns
        if i == 1
            fprintf(foutc,'if (p<cumpropensities[%d]) {\n',i-1);
        else
            fprintf(foutc,'} else if (p<cumpropensities[%d]) {\n',i-1);
        end
        fprintf(foutc,'  // rxn: %s\n',lines{i});
        allstoichiometry = rhsstoichiometry-lhsstoichiometry;
        idx = find(allstoichiometry(:,i)~=0);
        for j = idx'
            fprintf(foutc,'  %s=%s + %d;\n',uniquetok{j},uniquetok{j},allstoichiometry(j,i));
        end

        fprintf(foutc,'\n');

        for currrxn = 1:numrxns
            fprintf(foutc,'  //update propensity for %s\n',lines{currrxn});
            fprintf(foutc,'  propensities[%d] = %s',currrxn-1,ratenames{currrxn});

            if ismember(currrxn,1:n_species) %Production
                for ion = 1:2
                    if ion ==1
                        fprintf(foutc,'*%s*%s',on_ratenames{currrxn},uniquetok{2*currrxn+n_species});
                    else
                        fprintf(foutc,'+%s*%s',ratenames{currrxn},uniquetok{2*currrxn+n_species-1});
                    end
                end
            elseif ismember(currrxn,n_species+1:2*n_species) %Degradation
                fprintf(foutc,'*%s', uniquetok{currrxn-n_species});
            elseif ismember(currrxn,2*n_species+1:3*n_species) %Burst on
                 jnew = currrxn-(2*n_species);
                 countjparmat = 0;
                 for jparmat = 1:nspecies/3
                     if parmat(jparmat,jnew) == 1 && countjparmat == 0 %Look for all species for which the production is dependent on abundance of j
                         fprintf(foutc,'*(pow(%s,%s)/(pow(%s,%s)+pow(%s,%s))',uniquetok{jparmat},n_ratenames{jparmat},...
                             k_ratenames{jparmat},n_ratenames{jparmat},uniquetok{jparmat},n_ratenames{jparmat});
                         countjparmat = countjparmat +1;
                     elseif parmat(jparmat,jnew) == 1
                         fprintf(foutc,'+pow(%s,%s)/(pow(%s,%s)+pow(%s,%s))',uniquetok{jparmat},n_ratenames{jparmat},...
                             k_ratenames{jparmat},n_ratenames{jparmat},uniquetok{jparmat},n_ratenames{jparmat});
                         countjparmat = countjparmat +1;
                     end
                     if jparmat == nspecies/3 && countjparmat ~= 0
                         fprintf(foutc,')');
                     end
                 end
                 if countjparmat == 0
                     fprintf(foutc,'*0');
                 end
                 fprintf(foutc,'*%s',uniquetok{2*jnew+n_species-1});
                 fprintf(foutc,'+%s*%s',b_ratenames{jnew},uniquetok{2*jnew+n_species-1});
            elseif ismember(currrxn,3*n_species+1:4*n_species) %Burst off
                jnew = currrxn-(n_species*3);
                fprintf(foutc,'*%s',uniquetok{2*jnew+n_species});

            end

            fprintf(foutc,';\n');
        end
    end

    fprintf(foutc,'}\n');

    while 1  %Write out the rest of the file
        tline = fgetl(finc);
        if ~ischar(tline),break,end
        fprintf(foutc,'%s\n',tline);
    end

    fclose(finc);
    fclose(foutc);

    fprintf('Compiling %smex.c... ',outputfile);
    name=sprintf('mex %shistomex.c',outputfile);
    tic;
    eval(name)
    t = toc;
    fprintf('elapsed time = %g seconds\n',t);
end
%**********************************************%

%Generate the parameter file

fid = fopen([outputfile,'params.m'],'w');
fprintf(fid,'%% Parameter file\n\n\n');
fprintf(fid,'%% Simulation parameter values\n\n');
fprintf(fid,'%% Maximum number of Gillespie steps\n');
fprintf(fid,'maxgillespiesteps = %d;\n\n',maxgillespiesteps);
fprintf(fid,'%% Initial time\n');
fprintf(fid,'currT = 0;\n\n\n');

fprintf(fid,'%% Initial values for species\n\n');
for i = 1:nspecies
    fprintf(fid,'%s;\n',initiallines{i});
end

fprintf(fid,'\n\n%% EVERYTHING BELOW IS BOOKKEEPING; DO NOT ALTER!\n');
fprintf(fid,'\n\n%% Reaction rates\n\n');
for i = 1:numrxns
    fprintf(fid,'%% Rxn: %s\n',lines{i});
    fprintf(fid,'%s = %g;\n',ratenames{i},rate(i));
end

for i = 1:length(on_ratenames)
    fprintf(fid,'%s = %g;\n',on_ratenames{i},on_rate(i));
end

for i = 1:length(b_ratenames)
    fprintf(fid,'%s = %g;\n',b_ratenames{i},b_rate(i));
end

for i = 1:length(k_ratenames)
    fprintf(fid,'%s = %g;\n',k_ratenames{i},k_rate(i));
end

for i = 1:length(n_ratenames)
    fprintf(fid,'%s = %g;\n',n_ratenames{i},n_rate(i));
end

fprintf(fid,'numrxns = %d;\n',numrxns);
fprintf(fid,'nspecies = %d;\n',nspecies);

fprintf(fid,'species = zeros(nspecies,1);\n');
for i = 1:nspecies
    fprintf(fid,'species(%d) = %s;\n',i,uniquetok{i});
end

fprintf(fid,'rates = zeros(numrxns,1);\n');

for i = 1:numrxns
    fprintf(fid,'rates(%d) = %s;\n',i,ratenames{i});
end

for d = 1:length(on_ratenames)
    fprintf(fid,'rates(%d) = %s;\n',i+d,on_ratenames{d});
end

for l = 1:length(b_ratenames)
    fprintf(fid,'rates(%d) = %s;\n',i+l+d,b_ratenames{l});
end

for j = 1:length(k_ratenames)
    fprintf(fid,'rates(%d) = %s;\n',i+l+j+d,k_ratenames{j});
end

for k = 1:length(n_ratenames)
    fprintf(fid,'rates(%d) = %s;\n',i+l+j+k+d,n_ratenames{k});
end

fprintf(fid,'y0 = zeros(nspecies,1);\n');
for i = 1:nspecies
    fprintf(fid,'y0(%d) = %s;\n',i,uniquetok{i});
end

%Intialize the propensities
fprintf(fid,'%% Intialize the propensities...\n');
for currrxn = 1:numrxns
    fprintf(fid,'propensity(%d) = %s',currrxn,ratenames{currrxn});
    
    if ismember(currrxn,1:n_species) %Production
        for ion = 1:2
            if ion ==1
                fprintf(fid,'*%s*%s',on_ratenames{currrxn},uniquetok{2*currrxn+n_species});
            else
                fprintf(fid,'+%s*%s',ratenames{currrxn},uniquetok{2*currrxn+n_species-1});
            end
        end
    elseif ismember(currrxn,n_species+1:2*n_species) %Degradation
        fprintf(fid,'*%s', uniquetok{currrxn-n_species});
    elseif ismember(currrxn,2*n_species+1:3*n_species) %Burst on
        jnew = currrxn-(2*n_species);
        countjparmat = 0;
        for jparmat = 1:nspecies/3
            if parmat(jparmat,jnew) == 1 && countjparmat == 0 %look for all species for which the production is dependent on abundance of j
                fprintf(fid,'*((%s^%s)/(%s^%s+%s^%s)',uniquetok{jparmat},n_ratenames{jparmat},...
                    k_ratenames{jparmat},n_ratenames{jparmat},uniquetok{jparmat},n_ratenames{jparmat});
                countjparmat = countjparmat +1;
            elseif parmat(jparmat,jnew) == 1
                fprintf(fid,'+(%s^%s)/(%s^%s+%s^%s)',uniquetok{jparmat},n_ratenames{jparmat},...
                    k_ratenames{jparmat},n_ratenames{jparmat},uniquetok{jparmat},n_ratenames{jparmat});
                countjparmat = countjparmat +1;
            end
            if jparmat == nspecies/3 && countjparmat ~= 0
                fprintf(fid,')');
            end
        end
        if countjparmat == 0
            fprintf(fid,'*0');
        end
        fprintf(fid,'*%s',uniquetok{2*jnew+n_species-1});
        fprintf(fid,'+%s*%s',b_ratenames{jnew},uniquetok{2*jnew+n_species-1});
    elseif ismember(currrxn,3*n_species+1:4*n_species) %Burst off
        jnew = currrxn-(n_species*3);
        fprintf(fid,'*%s',uniquetok{2*jnew+n_species});      
    end
    
    fprintf(fid,';\n');
    
end

fprintf(fid,'\n');

fclose(fid);


end

