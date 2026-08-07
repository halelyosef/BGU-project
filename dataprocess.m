% List of input files - update paths as needed
%יש להוריד את העמדוה השנייה של שמות הרגיסטרים ואת החזרות בתדרים 
inputFiles = {
   
};

% Output folder - make sure it exists
outputFolder = ;
if ~isdir(outputFolder)
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
end

fprintf('=== Processing %d files to create cleaned versions ===\n\n', numel(inputFiles));

for fileIdx = 1:numel(inputFiles)
    inputFile = inputFiles{fileIdx};
    fprintf('[%d/%d] Processing: %s\n', fileIdx, numel(inputFiles), inputFile);
    
    try
        % ========== STEP 1: Read the file ==========
        rawData = readcell(inputFile);
        if isempty(rawData)
            warning('  File is empty. Skipping.');
            continue;
        end
        fprintf('  Raw size: [%d x %d]\n', size(rawData,1), size(rawData,2));
        
        % ========== STEP 2: Process header row (row 1, columns 2:end) ==========
        headerRaw = rawData(1,2:end);
        % Convert to string, remove 'freq', keep only digits and dot
        headerStr = string(headerRaw);
        headerStr = erase(headerStr, "freq");
        headerStr = regexprep(headerStr, '[^\d.]', ''); % keep only digits and dot
        headerNum = str2double(headerStr);
        % Check for conversion errors
        invalidNum = isnan(headerNum);
        if any(invalidNum)
            fprintf('  Warning: %d header values could not be converted to numbers.\n', sum(invalidNum));
            headerNum(invalidNum) = 0; % replace invalid with 0
        end
        headerMHz = headerNum ./ 1e6; % convert to MHz
        % Format as clean string (avoid scientific notation, trim)
        headerStrClean = string(arrayfun(@(x) compose('%.6g', x), headerMHz, 'UniformOutput', false));
        headerStrClean = strtrim(headerStrClean);
        % Update rawData
        data = rawData;
        data(1,2:end) = num2cell(headerStrClean);
        
        % ========== STEP 3: Process first column: remove '0x' and trim ==========
        addrRaw = rawData(2:end,1);
        addrStr = string(addrRaw);
        addrStr = erase(addrStr, "0x");
        addrStr = strtrim(addrStr);
        % Replace empty strings with placeholder
        addrStr(addrStr == '') = '<EMPTY_ADDR>';
        data(2:end,1) = num2cell(addrStr);
        
        % ========== STEP 4: Replace missing values with 0 ==========
        data(ismissing(data)) = {0};
        
        % ========== STEP 5: Build output filename ==========
        [~, name, ext] = fileparts(inputFile);
        % Remove any problematic chars from name for safety
        safeName = regexprep(name, '[\\\/\:\*\?\"\<\>\|]', '_');
        outputFile = fullfile(outputFolder, [safeName, '_cleaned.xlsx']); % ✅ FIXED: fullfile, not fullpath
        
        % ========== STEP 6: Save cleaned data ==========
        writecell(data, outputFile);
        fprintf('  ✅ Saved cleaned file to: %s\n', outputFile);
        fprintf('     Size: [%d rows × %d columns]\n', size(data,1), size(data,2));
        
    catch ME
        fprintf('  ❌ ERROR processing %s: %s\n', inputFile, ME.message);
        fprintf('     Skipping this file.\n');
    end
    
    fprintf('\n'); % blank line for readability
end

fprintf('🎉 All files processed. Cleaned files are in: %s\n', outputFolder);


inputFiles = {
   
};

outputFolder = 'L:\6412\Users\project maor and halel\MEAS\CLEAN\1MHz';

outputFile = "merged_clean.xlsx";

if ~exist(outputFolder,'dir');
    mkdir(outputFolder);
end

fullOutputPath = fullfile(outputFolder, outputFile);


%% ===== READ FIRST FILE (BASE) =====
Tbase = readtable(inputFiles{1}, 'PreserveVariableNames', true);

% first column is registers
regBase = Tbase{:,1};
freqVarsBase = Tbase.Properties.VariableNames(2:end);

%% ===== LOOP OVER OTHER FILES =====
for f = 2:numel(inputFiles)

    T = readtable(inputFiles{f}, 'PreserveVariableNames', true);

    reg = T{:,1};
    freqVars = T.Properties.VariableNames(2:end);

    % check register consistency
    if ~isequal(reg, regBase)
        error('Register list mismatch in file: %s', inputFiles{f});
    end

    % identify new frequencies (columns)
    newFreqVars = setdiff(freqVars, freqVarsBase, 'stable');

    if isempty(newFreqVars)
        continue; % nothing to add
    end

    % append only new frequency columns
    Tbase = [Tbase, T(:, newFreqVars)];

    % update base frequency list
    freqVarsBase = Tbase.Properties.VariableNames(2:end);
end

%% ===== SORT FREQUENCIES NUMERICALLY =====
% extract numeric freq values
freqNums = cellfun(@str2double, freqVarsBase);

[~, idx] = sort(freqNums);

Tbase = [Tbase(:,1), Tbase(:, idx+1)];

%% ===== WRITE OUTPUT =====
writetable(Tbase, fullOutputPath);
disp(['Merged file saved to:', fullOutputPath]);

%% ================= USER SETTINGS =================

inputFile  = ';
outputFile = ;


bitTable = {
    '0x271', 7, 0
    '0x272', 2, 0
    '0x273', 7, 0
    '0x274', 7, 0
    '0x275', 6, 0
    '0x282', 2, 0
    '0x279', 3, 0
    '0x282', 4, 3
    '0x27B', 5, 0
    '0x280', 3, 0
    '0x27F', 3, 0
    '0x27E', 3, 0
    '0x27E', 7, 4
    '0x27F', 7, 4
    '0x290', 6, 4
    '0x005', 7, 4
    '0x278', 6, 3
    '0x291', 3, 0
    '0x277', 7, 0
    '0x276', 6, 0
    '0x278', 0, 0
};
%% ================= READ EXCEL ====================
T = readtable(inputFile, 'PreserveVariableNames', true);

% Extract frequencies (row 1, columns 2:end)
freqs = T.Properties.VariableNames(2:end);

% Register names (column 1, rows 1:end)
regs = string(T{:,1});

% Numeric data
data = T{:,2:end};

%% ================= PREP BIT TABLE ================
bitRegs = string(bitTable(:,1));
bitRegs = erase(bitRegs, "0x");          % remove 0x
bitRegs = upper(bitRegs);

msbList = cell2mat(bitTable(:,2));
lsbList = cell2mat(bitTable(:,3));

%% ================= PROCESS =======================
outRegs  = strings(size(bitTable,1),1);
outData  = nan(size(bitTable,1), size(data,2));

for k = 1:size(bitTable,1)

    regName = bitRegs(k);
    msb = msbList(k);
    lsb = lsbList(k);

    % Find register row
    idx = find(upper(regs) == regName, 1);

    if isempty(idx)
        warning("Register %s not found in input", regName);
        continue
    end

    % Output register label
    if msb == lsb
        outRegs(k) = regName + "[" + msb + "]";
    else
        outRegs(k) = regName + "[" + msb + ":" + lsb + "]";
    end

    for c = 1:size(data,2)
        val = data(idx,c);

        if isnan(val)
            continue
        end

        % Decimal → binary (8 bits)
        bin = dec2bin(val, 8) - '0';   % vector [b7 ... b0]

        % Extract bits
        slice = bin(8-msb : 8-lsb);

        % Binary → decimal
        outData(k,c) = bin2dec(char(slice + '0'));
    end
end

%% ================= WRITE OUTPUT ==================
outTable = array2table(outData, 'VariableNames', freqs);
outTable = addvars(outTable, outRegs, 'Before', 1, ...
                   'NewVariableNames', 'reg_bits');

writetable(outTable, outputFile);

% %%=======PLOTTING============
% freqHz = zeros(1, numel(freqs));
% for i = 1:numel(freqs)
%     freqHz(i) = str2double(freqs{i});
% end 
% 
% for r= 1:size(outData,1)
%     y = outData(r, :);
%     if all(isnan(y));
%         continue
%     end 
% 
%     figure;
%     plot(freqHz, y, 'o', 'LineWidth',1.5);
%     grid on; 
% 
%     xlabel('freq [MHz]');
%     ylabel('reg');
%     title(outRegs(r), 'Interpreter','none');
% end 