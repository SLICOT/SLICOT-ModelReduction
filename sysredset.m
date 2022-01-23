function options = sysredset(varargin)
%SYSREDSET Create/alter SYSRED OPTIONS structure.
%   OPTIONS = SYSREDSET('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a
%   system order reduction options structure OPTIONS in which the named 
%   parameters have the specified values.  Any unspecified parameters
%   are set to [] (a parameter with value [] indicate to use the default value 
%   for that parameter when OPTIONS is passed to the system reduction function). 
%   It is sufficient to type only the leading characters that uniquely identify 
%   the parameter.  Case is ignored for parameter names.  
%   NOTE: For values that are strings, correct case and the complete string  
%   are required; if an invalid string is provided, the default is used.
%   
%   OPTIONS = SYSREDSET(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS 
%   with the named parameters altered with the specified values.
%   
%   OPTIONS = SYSREDSET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in 
%   OLDOPTS. 
%   
%   SYSREDSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values, with defaults shown in {} 
%   when the default is the same for all functions that use that option.
%   Use SYSREDSET(SYSREDFUNCTION) to see options for a specific function.
%
%   OPTIONS = SYSREDSET (with no input arguments) creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   OPTIONS = SYSREDSET(SYSREDFUNCTION) creates an options structure with all
%   the parameter names and default values relevant to the system reduction
%   function named in SYSREDFUNCTION. For example,
%           SYSREDSET('bst') 
%   or
%           SYSREDSET(@bst)
%   returns an options structure containing all the parameter names and  
%   default values relevant to the function 'bst'.
%   
%SYSREDSET PARAMETERS
%BalredMethod - Balancing reduction approach used: the Balanced Truncation Approximation
%             (BTA) or the Singular Perturbation Approximation (SPA) [ {bta} | spa ]
%AccuracyEnhancing - Accuracy enhancing technique used: the balancing-free square-root
%             (BFSR) method or square-root (SR) method [ {bfsr} | sr ]
%Tolred     - Tolerance on Hankel-singular values to determine the order of the
%             reduced order model [ positive scalar {0} ]; if Tolred = 0, a value of Tolred
%             appropriate for mininimal realization is internally determined
%TolMinreal - Tolerance on Hankel-singular values to determine the order of a minimal
%             realization [ positive scalar {0} ]; if TolMinreal = 0, a value of TolMinreal
%             is internally determined
%CStabDeg   - Stability degree parameter for continuous-time systems
%             [ nonpositive scalar | {-sqrt(eps)} ]
%DStabDeg   - Stability degree parameter for discrete-time systems
%             [ nonpositive scalar | {1-sqrt(eps)} ]
%Order      - Order of the reduced system [ integer {-1} ]; Order = -1 means that the order
%             is determined automatically using the value of Tolred.
%BstBeta    - Weighting factor for the BST method [ scalar {0}  ]; BstBETA = 0 means
%             pure relative error method (BST), while BstBETA > 0 means weighted 
%             relative-absolute (BST-BTA) method.
%FWEContrGramian - Choice of controllability Grammian for the frequency-weighted balancing 
%             related method [ {standard} | enhanced ]
%FWEObservGramian - Choice of observability Grammian for the frequency-weighted balancing 
%             related method [ {standard} | enhanced ]
%FWEAlphaContr - Weighting factor for the frequency-weighted controllability Grammian 
%             [ positive subunitary scalar {0} ]
%FWEAlphaObserv - Weighting factor for the frequency-weighted observability Grammian 
%             [ positive subunitary scalar {0} ]
%FWEConredMethod - Choice of weighting factors for the frequency-weighted controller
%             reduction [ none | outputstab | inputstab | {performance} ]
%FWEHNAMethod - Choice of method to handle system inverses
%             [ {auto}  | inverse | noinverse ]
%FWEHNAopV  - Operation to be performed on left weight V [ {none} | inv | conj | cinv ] 
%FWEHNAopW  - Operation to be performed on right weight W [ {none} | inv | conj | cinv ] 
%FWEOptimize  - Optimization to be performed on computed reduced model [ none | {d} | cd ] 
%CFConredMethod - Coprime factorization based controller reduction method [ {fwe}  | nofwe ]
%CoprimeFactorization - Type of coprime factorization [ left | {right} ]
%
%   See also SYSREDGET.

%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('               Function: [ {balred} | hna | bst | cof  | fwbred | fwhna | fwbconred | cfconred ]\n');
    fprintf('           BalredMethod: [ {bta} | spa ]\n');
    fprintf('      AccuracyEnhancing: [ {bfsr} | sr ]\n');
    fprintf('                 Tolred: [ positive scalar {0} ]\n');
    fprintf('             TolMinreal: [ positive scalar {0} ]\n');
    fprintf('                  Order: [ integer {-1} ]\n');
    fprintf('               CStabDeg: [ nonpositive scalar | {-sqrt(eps)} ]\n');
    fprintf('               DStabDeg: [ subunitary scalar | {1-sqrt(eps)} ]\n');
    fprintf('                BstBeta: [ scalar {0}  ]\n');
    fprintf('       FWEContrGramian: [ {standard} | enhanced ]\n');
    fprintf('      FWEObservGramian: [ {standard} | enhanced ]\n');
    fprintf('          FWEAlphaContr: [ positive subunitary scalar {0} ]\n');
    fprintf('         FWEAlphaObserv: [ positive subunitary scalar {0} ]\n');
    fprintf('   CoprimeFactorization: [ left | {right} ]\n');
    fprintf('           OutputWeight: [ {stab} | perf | none]\n');
    fprintf('            InputWeight: [ {stab} | none]\n');
    fprintf('         CFConredMethod: [ {fwe}  | nofwe ]\n');
    fprintf('        FWEConredMethod: [ none | outputstab | inputstab | {performance} ]\n');
    fprintf('           FWEHNAMethod: [ {auto}  | inverse | noinverse ]\n');
    fprintf('              FWEHNAopV: [ {none} | inv | conj | cinv ]\n'); 
    fprintf('              FWEHNAopW: [ {none} | inv | conj | cinv ]\n'); 
    fprintf('            FWEOptimize: [ none | {d} | cd ]\n'); 
    fprintf('\n');
    return;
end

options = struct(  ...
    'BalredMethod', [], ...
    'AccuracyEnhancing', [], ...
    'TolRed', [], ...
    'TolMinreal', [], ...
    'Order', [], ...
    'CStabDeg', [], ...
    'DStabDeg', [], ...
    'BstBeta', [], ...
    'FWEContrGramian', [], ...
    'FWEObservGramian', [], ...
    'FWEAlphaContr', [], ...
    'FWEAlphaObserv', [], ...
    'CoprimeFactorization', [], ...
    'OutputWeight', [], ...
    'InputWeight', [], ...
    'CFConredMethod', [], ...
    'FWEConredMethod', [], ...
    'FWEHNAMethod', [], ...
    'FWEHNAopV', [], ...
    'FWEHNAopW', [], ...
    'FWEOptimize', []);

numberargs = nargin; % we might change this value, so assign it

% If we pass in a function name then return the defaults.
if (numberargs == 1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            error('No default options available: the function ''%s'' does not exist on the path.',funcname)
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch
        error('No default options available for the function ''%s''.',funcname)
    end
    % To get output, run the rest of SYSREDSET as if called with SYSREDSET(options, optionsfcn)
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = fieldnames(options);
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(['Expected argument %d to be a string parameter name ' ...
                    'or an options structure\ncreated with SYSREDSET.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = getfield(arg, Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                [valid, errmsg] = checkfield(Names{j,:},val);
                if valid
                    options = setfield(options, Names{j,:},val);
                else
                    error(errmsg);
                end
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error('Expected argument %d to be a string parameter name.', i);
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error('Unrecognized parameter name ''%s''.', arg);
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                error('%s).', msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else           
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        [valid, errmsg] = checkfield(Names{j,:},arg);
        if valid
            options = setfield(options, Names{j,:},arg);
        else
            error(errmsg);
        end
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error('Expected value for parameter ''%s''.', arg);
end

%-------------------------------------------------
function f = getfield(s,field)
%GETFIELD Get structure field contents.
%   F = GETFIELD(S,'field') returns the contents of the specified
%   field.  This is equivalent to the syntax F = S.field.
%   S must be a 1-by-1 structure.  
% 

sref.type = '.'; sref.subs = field;
f = subsref(s,sref);

%-------------------------------------------------
function s = setfield(s,field,value)
%SETFIELD Set structure field contents.
%   S = SETFIELD(S,'field',V) sets the contents of the specified
%   field to the value V.  This is equivalent to the syntax S.field = V.
%   S must be a 1-by-1 structure.  The changed structure is returned.
%

sref.type = '.'; sref.subs = field;
s = subsasgn(s,sref,value);

%-------------------------------------------------
function [valid, errmsg] = checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   [VALID, MSG] = CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%


valid = 1;
errmsg = '';
% empty matrix is always valid
if isempty(value)
    return
end

switch field
case {'TolRed','TolMinreal','BstBeta'} % real positive scalar
    if ~(isa(value,'double') && value >= 0),
        valid = 0;
        if ischar(value)
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive number (not a string).',field);
        else
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive number.',field);
        end
    end
case {'Order'} % real positive scalar or -1
    if ~(isa(value,'double') && (value >= 0 || value == -1 ) ),
        valid = 0;
        if ischar(value)
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive number (not a string).',field);
        else
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive number.',field);
        end
    end
case {'CStabDeg'} % real negative scalar
    if ~(isa(value,'double') && value <= 0) 
        valid = 0;
        if ischar(value)
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real nonpositive number (not a string).',field);
        else
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real nonpositive number.',field);
        end
    end
case {'DStabDeg','FWEAlphaContr','FWEAlphaObserv'} % real subunitary positive scalar
    if ~(isa(value,'double') && value <= 1 && value >= 0) 
        valid = 0;
        if ischar(value)
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive subunitary number (not a string).',field);
        else
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive subunitary number.',field);
        end
    end
case {'BalredMethod'} % bta, spa
    if ~isa(value,'char') || ~any(strcmp(value,{'bta';'spa';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''bta'', or ''spa''.',field);
    end
case {'AccuracyEnhancing'} % bfsr, sr
    if ~isa(value,'char') || ~any(strcmp(value,{'bfsr';'sr';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''bfsr'', or ''sr''.',field);
    end
case {'FWEContrGramian','FWEObservGramian'} % standard, enhanced 
    if ~isa(value,'char') || ~any(strcmp(value,{'standard';'enhanced';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''standard'', or ''enhanced''.',field);
    end
case {'FWEConredMethod'} % none, outputstab, inputstab, performance
    if ~isa(value,'char') || ~any(strcmp(value,{'none';'outputstab';'inputstab';'performance';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'', ''outputstab'' , ''inputstab'' or ''performance''.',field);
    end
case {'FWEHNAMethod'} % auto, inverse, noinverse 
    if ~isa(value,'char') || ~any(strcmp(value,{'auto';'inverse';'noinverse';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''auto'', ''inverse'' or ''noinverse''.',field);
    end
case {'FWEHNAopV','FWEHNAopW'} % none, inv, conj, cinv  
    if ~isa(value,'char') || ~any(strcmp(value,{'none';'inv';'conj';'cinv';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'', ''inv'', ''conj'' or ''cinv''.',field);
    end
case {'FWEOptimize'} % none, b, c, d  
    if ~isa(value,'char') || ~any(strcmp(value,{'none';'cd';'d';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'', ''b'', ''bd'', ''c'', ''cd'' or ''d''.',field);
    end
case {'CFConredMethod'} % fwe, nofwe 
    if ~isa(value,'char') || ~any(strcmp(value,{'fwe';'nofwe';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''fwe'' or ''nofwe''.',field);
    end
case {'CoprimeFactorization'} % left, right 
    if ~isa(value,'char') || ~any(strcmp(value,{'left';'right';}))
        valid = 0;
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''left'' or ''right''.',field);
    end
otherwise  
    valid = 0;
    error('Unknown field name for Options structure.')
end
