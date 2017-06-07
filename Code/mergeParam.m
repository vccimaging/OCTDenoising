function newpar = mergeParam(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEWPARAM = MERGEPARAM(PARAM1,PARAM2,...)
% MERGEPARAM merges parameter sets PARAM1, PARAM2,... into NEWPARAM.
% Each parameter set as well as the new one is a struct with field names. 
% If a field is defined multiple times, the last input argument is 
% taken into account.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newpar = varargin{1};
if ~isstruct(newpar);
    error('mergeParam:1st parameter set is not a structure');
end

for i = 2:length(varargin)
    temp = varargin{i};
    if ~isempty(temp)
        if ~isstruct(temp);
            
            switch i
                case 2
                    error('mergeParam:2nd parameter set is not a structure');
                case 3
                    error('mergeParam:3rd parameter set is not a structure');
                otherwise
                    error(['mergeParam:' num2str(i) '-th parameter set is not a structure']);
            end
        end
        for f = fieldnames(temp)'
            %newparam = setfield(newparam,f{1},getfield(temp,f{1}));
            newpar.(f{1}) = temp.(f{1});
        end
    end
end

