function [] = main_configstepping(inputfile,outputfile)
%MAIN_CONFIGSTEPPING Assemble resistance matrices for a sequence of configurations.
%
%   MAIN_CONFIGSTEPPING(inputfile,outputfile) loads a set of filament
%   configurations from INPUTFILE, computes the full slender-body theory
%   interaction matrix for each configuration, converts it to a resistance
%   matrix, and saves the results to OUTPUTFILE.
%
%   This routine supports restarting from partial output when the input
%   file status is marked as 'inprogress' and OUTPUTFILE already exists.
%
%   INPUTS
%       inputfile   MAT-file containing simulation/configuration data.
%                   Required variables include:
%                       status, epsil, psi, Nturns, c, NLegendre,
%                       rtolr, atolr, x, es
%       outputfile  MAT-file used to save intermediate and final results
%
%   OUTPUTS
%       No direct output arguments.
%       Results are written to OUTPUTFILE, including:
%           FullRes   Full resistance matrix for each configuration
%           ResM      Single-filament resistance matrix
%           SbtM      Single-filament SBT matrix
%           fmom      Force-moment quantity returned by SBTSELF
%
%   NOTES
%       For each configuration, the routine:
%           1. Rotates the single-filament SBT matrix into each filament frame
%           2. Builds the block interaction matrix including pairwise terms
%           3. Converts the full SBT matrix to a resistance matrix
%
%       The third dimension of x and es is assumed to index filaments, and
%       the fourth dimension is assumed to index configurations/time steps.
%
%   EXAMPLE
%       main_configstepping('input_configs.mat','output_results.mat')

% Load geometric, numerical, and configuration data from the input file.
disp(['loading ' inputfile])
load(inputfile,'status','epsil','psi','Nturns','c','NLegendre','rtolr','atolr',...
    'x','es')

% Number of filaments is stored in the third dimension.
Nfilaments = size(x,3); 

% Load previously computed single-filament quantities if available.
% If a run is being resumed and the output file already exists, prefer the
% saved copies there; otherwise load them from the input file.
if (strcmp(status,'inprogress') && isfile(outputfile))
    load(outputfile,'SbtM')
else
    load(inputfile,'SbtM')
end

% Compute the isolated-filament resistance/SBT data only if it has not
% already been stored.
if ~exist('SbtM','var')
    [~, SbtM, ~] = sbtself(epsil,psi,Nturns,c,NLegendre,rtolr,atolr);
end

% Decide whether to resume an existing run or start from scratch.
switch true
    case (strcmp(status,'inprogress') && isfile(outputfile))
        % Resume from the last saved configuration index.
        load(outputfile,'kk','FullRes')
        start_kk = kk;

    case strcmp(status,'notstarted')
        % Initialise storage for the full resistance matrix at each
        % configuration and mark the input file as in progress.
        FullRes = zeros(Nfilaments*6,Nfilaments*6,size(x,4));
        start_kk = 1;
        status = 'inprogress';
        save(inputfile,'status','-append')
end

% Start timing the current run.
starttime = cputime;

% Loop over all stored configurations.
for kk = start_kk:size(x,4)
    % Report progress and save intermediate state for restart capability.
    disp([inputfile ' progress:' num2str(kk) '/' num2str(size(x,4)) ...
        ', running time: ' num2str((cputime-starttime)/3600) 'h'])
    save(outputfile)
    
    % Extract positions and orientation frames for the current configuration.
    x0 = x(:,:,:,kk);
    es0 = es(:,:,:,kk);
    
    % Rotate the single-filament SBT matrix into the lab-frame orientation
    % of each filament.
    SBTnn = zeros([size(SbtM) Nfilaments]);
    for n=1:Nfilaments
       SBTnn(:,:,n) = rotateSBT(SbtM,es0(:,:,n));
    end
    
    % Initialise the full block SBT matrix for this configuration.
    FullSBT = zeros(Nfilaments*size(SBTnn,1),Nfilaments*size(SBTnn,2));
    
    % Fill diagonal blocks with self-interaction matrices.
    for n=1:Nfilaments
       rowloc = size(SBTnn,1)*(n-1)+(1:size(SBTnn,1));
       colloc = size(SBTnn,2)*(n-1)+(1:size(SBTnn,2));
       FullSBT(rowloc,colloc) = SBTnn(:,:,n);
    end
    
    % Fill off-diagonal blocks with pairwise hydrodynamic interactions.
    % Only compute each pair once, then use symmetry to fill the transpose.
    for n=1:Nfilaments
        for m=(n+1):Nfilaments
            % Compute cross-interaction matrix
            SBTnm = sbtcross(x0(:,:,n),x0(:,:,m),es0(:,:,n),es0(:,:,m),...
                epsil,psi,Nturns,c,NLegendre,rtolr,atolr);
            rowloc = size(SBTnm,1)*(n-1)+(1:size(SBTnm,1)); % nth row (of block matrices)
            colloc = size(SBTnm,2)*(m-1)+(1:size(SBTnm,2)); % mth col (of block matrices)
            FullSBT(rowloc,colloc) = SBTnm;
            
            % Reverse interaction block obtained by transpose symmetry.
            FullSBT(colloc,rowloc) = SBTnm';
        end
    end
    
    % Convert the full SBT matrix into the resistance matrix for this
    % configuration.
    FullRes(:,:,kk) = sbt2res(FullSBT,Nfilaments,es0,psi,Nturns,c,NLegendre);
end

% Record total CPU time for this run.
runtime = cputime-starttime;

% Save final results to the output file.
save(outputfile)
disp(['saving ' outputfile ', total run time: ' num2str(runtime/3600) 'h'])

% Mark the input file as fully processed.
status = 'done';
save(inputfile,'status','-append')
end