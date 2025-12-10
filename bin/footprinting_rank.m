%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Footrpints Detector v1.0
% Copyright (c) 2008 by Xiaoyu Chen and William Noble
% All rights reserved.
% Redistribution is not permitted without the permission of the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function footprinting_rank(DHSCutCountsFile, ...
                         mappableFile, ...
                         minK, ...
                         maxK, ...
                         KStep, ...
                         isDHS, ...
                         localBgrWidth, ...
                         flankSize, ...
                         maxColNum, ...
                         rankThresh, ...
                         qValueThresh, ...
                         transformType, ...
                         isOneSideBkg, ...
                         outPath, ...
                         outPrefix)

% Check if statistics package is loaded
if ~exist('normcdf')
    error('Statistics package is required. Please install it using: pkg install -forge statistics');
end

% localBgrWidth must be odd!
if flankSize < floor(localBgrWidth/2)
  error('flank-size must be larger than half of local-background-width!');
end

if outPath(end) ~= '/'
  outPath = [outPath '/'];
end

fpCC = fopen(DHSCutCountsFile, 'r');
fpMap = fopen(mappableFile, 'r');
fpDataTemp = fopen([outPath outPrefix '.datatemp'], 'w');
fpNullTemp = fopen([outPath outPrefix '.nulltemp'], 'w');

% reset seed using rng for compatibility
rng(99999999);
rowCount = 0;

if localBgrWidth ~= 0
  [winStart, winStop] = window_positions(minK, maxK, KStep, maxColNum);
end

while not(feof(fpCC))
  rowCount = rowCount + 1; 
  if mod(rowCount, 1000) == 1
    fprintf('row %d \n', rowCount);
  end

  line = fgetl(fpCC);
  DHSCC = str2num(line);
  ix = find(isnan(DHSCC));
  DHSCC(ix) = 0;

  line = fgetl(fpMap);
  mappable = str2num(line);
 
  % This fillter ONLY works for human DHSs!
  if (isDHS) && (mean(DHSCC((1+flankSize):(end-flankSize))) <= 1)
    continue;
  end 

  ix = find(mappable ~= 1); 
  mappable(ix) = 0;


  %% ONLY for human DHSs, filter ending regions with low cut-counts
  %if isDHS
  %  [realStart, realStop] = filter_ends(DHSCC, flankSize);
  %else
  
  realStart = 1 + flankSize;
  realStop = min(length(DHSCC) - flankSize, maxColNum);
  %end

  % now compute pvalues for each kmer
  compute_pvalues_local(rowCount, DHSCC, mappable,...
                        realStart, realStop, transformType, ...
                        minK, maxK, KStep, localBgrWidth, ...
                        flankSize, winStart, winStop, ...
                        isOneSideBkg, ...
                        fpDataTemp, fpNullTemp);
end

fclose(fpMap);
fclose(fpCC);
fclose(fpDataTemp);
fclose(fpNullTemp);



tmp = importdata([outPath outPrefix '.datatemp']);
finalDHSId = tmp(:, 1);
finalStart = tmp(:, 2);
finalStop = tmp(:, 3);
finalPV = tmp(:, 4);
fprintf('%d candidates of footprints identified from data.\n', length(finalPV));

tmp = importdata([outPath outPrefix '.nulltemp']);
nullPV = tmp(:, 2);
clear tmp;
fprintf('%d candidates of footprints identified from null.\n', length(nullPV));

% compute qvalues directly based on nullpv
finalQV = compute_QV(finalPV, nullPV);

% output
fprintf('output to files...\n');
fpOutQV = fopen([outPath outPrefix '.qv' num2str(qValueThresh)], 'w');
ix = find(finalQV <= qValueThresh);
fprintf('%d significant footprints identified.\n', length(ix));
for i = 1:length(ix)
  id = finalDHSId(ix(i));
  start = finalStart(ix(i));
  stop = finalStop(ix(i));
  fprintf(fpOutQV, '%d\t%d\t%d\t%g\t%g\n', ...
          id, start, stop, finalPV(ix(i)), finalQV(ix(i)));
end
fclose(fpOutQV);



%function 
% filter regions that are right besides the ends and with low cut-counts
function [DHSStart, DHSStop] = filter_ends(data, flankSize)
firstCol = 1 + flankSize;
lastCol = length(data) - flankSize;

sumThresh = 0.05*sum(data(firstCol:lastCol));    
filter = 1;

currPos = firstCol;
while (filter == 1) && (currPos < lastCol)
  currSum = sum(data(firstCol:currPos)); 
  if currSum > sumThresh
    filter = 0;
  end
  currPos = currPos + 1;
end
DHSStart = currPos;
    

filter = 1;
currPos = lastCol;
while (filter == 1) && (currPos > DHSStart)
  currSum = sum(data(currPos:lastCol)); 
  if currSum > sumThresh
    filter = 0;
  end
  currPos = currPos - 1;
end
DHSStop = currPos;



% function
function [inputData] = generate_input(DHSCC, ...
                                      mappable, ...
                                      realStart, ...
                                      realStop, ...
                                      transformType)

cnum = length(DHSCC);
% converting cut-counts to ranks
if transformType == 1
  cutRanks = zeros(1, cnum);
  rk1 = [];
  rk2 = [];
 
  validIx = find(mappable(realStart:realStop) ~= 0) + realStart - 1;   
  if length(validIx) > 0
    validItem = DHSCC(validIx) ;
    [junk, ix1] = sort(validItem);
    rk1(ix1) = 1:length(ix1);
    [junk, ix2] = sort(validItem(end:-1:1));
    rk2(ix2) = 1:length(ix2);
    cutRanks(validIx) = (rk1 + rk2(end:-1:1))/2;
  end
  inputData = cutRanks;    

% truncate top 2%
elseif transformType == 2
  inputData = zeros(1, cnum);;
  validIx = find(mappable(realStart:realStop) ~= 0) + realStart -1;  
  if length(validIx) > 0
    validItem = DHSCC(validIx);
    currPrt = prctile(validItem, 95);
    validItem(find(validItem > currPrt)) = floor(currPrt);
    inputData(validIx) = validItem;
  end

% log(.+1) 
elseif transformType == 3
  inputData = ceil(log(DHSCC+1));
end



% function
% Remove kmers that overlapping with those of higher pvalues
function [cleanSc, cleanStart, cleanStop] = remove_overlaps(...
                                            score, ...
                                            start, ...
                                            stop)
flag = repmat(1, size(score));
[sortedSc, sortedIx] = sort(score);
sortedStart = start(sortedIx);
sortedStop = stop(sortedIx);

for i = 1:length(sortedSc)
  if flag(i) == 1
    for j = i+1:length(sortedSc)
      if (sortedStart(j) <= sortedStop(i) && sortedStop(j) >= sortedStart(i))...
         || (sortedStart(i) <= sortedStop(j) && sortedStop(i) >= sortedStart(j))
        flag(j) = 0;
      end
    end
  end
end

goodIx = find(flag == 1);
cleanSc = sortedSc(goodIx);
cleanStart = sortedStart(goodIx);
cleanStop = sortedStop(goodIx);



% function
% compute qvalues (fdr) from nullpv
function finalQV = compute_QV(finalSc, nullSc)

f = length(finalSc) / length(nullSc);

finalQV = [];
[sortedSc, sortedIx] = sort(finalSc);
numSamples = length(finalSc);

falsePos = length(find(nullSc <= sortedSc(numSamples)));
finalQV(sortedIx(numSamples)) = (f * falsePos) / numSamples;
lastSc = sortedSc(numSamples);

for j = (numSamples-1):-1:1
  if lastSc == sortedSc(j)
     finalQV(sortedIx(j)) = finalQV(sortedIx(j+1));
  else
     falsePos = length(find(nullSc <= sortedSc(j)));
     fdr = (f * falsePos) / j;
     finalQV(sortedIx(j)) = min(fdr, finalQV(sortedIx(j+1)));
     lastSc = sortedSc(j);
  end
end



% function
function [winStart, winStop] = window_positions(minK, maxK, KStep, maxCol)

winStart = zeros(maxK-minK+1, maxCol);
winStop = zeros(maxK-minK+1, maxCol);

for k = minK:KStep:maxK
  right = floor(k / 2);
  left = k - right -1;
  for c = 1:maxCol
    winStart(k,c) = c - left;
    winStop(k,c) = c + right;
  end
end



% function
function [kmerPV] = compute_PV(p, kmerCC, bkgCC)            

if (bkgCC*p > 10) && (bkgCC*(1-p) > 10)
  u = bkgCC * p;
  sigm = sqrt(bkgCC * p * (1-p));
  kmerPV = normcdf(kmerCC, u, sigm);
else
  kmerPV = binocdf(kmerCC, bkgCC, p);
end



% function
function compute_pvalues_local(rowNumber, data, mappable, ...
                               realStart, realStop, transformType, ...
                               minK, maxK, KStep, localBkWidth, ...
                               flankSize, winStart, winStop, ...
                               isOneSideBkg, ...
                               fpDataTemp, fpNullTemp)


% bkLen is required to be odd.
bkLen = localBkWidth;
bkHalf = floor(bkLen / 2);
firstCol = 1 + flankSize;
lastCol = length(data) - flankSize;
numK = floor((maxK - minK) / KStep) + 1;

mappableIx = find(mappable(realStart:realStop) ~= 0) + realStart - 1;
validLen = length(mappableIx);

if (validLen > minK)
  kmMaxNum = (realStop - realStart + 1) * numK;
  allPV = repmat(-1, 1, kmMaxNum);
  allStart = repmat(-1, 1, kmMaxNum);
  allStop = repmat(-1, 1, kmMaxNum);
  kmAllCount = 1;

  for c = realStart:realStop
    currBkStart = c - bkHalf;
    currBkStop = c + bkHalf;
   
    currBkMappable = mappable(currBkStart:currBkStop);
    if transformType == 0
      currBk = data(currBkStart:currBkStop);
    else
      currBk = generate_input(data(currBkStart:currBkStop), ...
                              currBkMappable,...
                              1, bkLen, transformType);
    end

    if ~isOneSideBkg
      currBkMapIx = find(currBkMappable ~= 0);
      bkRealW = sum(currBkMappable); 
      bkRealCC = floor(sum(currBk(currBkMapIx)));
    end

    % different length of kmer
    for k = minK:KStep:maxK
   
      if (winStart(k,c) >= firstCol) && (winStop(k,c) <= lastCol) 
        localStart = winStart(k,c) - currBkStart + 1;
        localStop = winStop(k,c) - currBkStart + 1;

        currMapIx = find(currBkMappable(localStart:localStop) ~= 0); 
        validK = length(currMapIx);

        if currBkMappable(localStart) && currBkMappable(localStop) && (validK >= minK)
          realK = sum(currBkMappable(localStart:localStop));
          currKmer = currBk(localStart:localStop);
          kmerCC = sum(currKmer(currMapIx));

          if ~isOneSideBkg
            p = realK/bkRealW;
            kmerPV = compute_PV(p, kmerCC, bkRealCC);
            allPV(kmAllCount) = kmerPV;
          else
            % left side background
            bkLeftW = sum(currBkMappable(1:localStop));
            currBkMapIx = find(currBkMappable(1:localStop) ~= 0);
            bkLeft = currBk(1:localStop);
            bkLeftCC = sum(bkLeft(currBkMapIx));

            p = realK/bkLeftW;
            kmerLeftPV = compute_PV(p, kmerCC, bkLeftCC);

            % right side background
            bkRightW = sum(currBkMappable(localStart:end));
            currBkMapIx = find(currBkMappable(localStart:end) ~= 0);
            bkRight = currBk(localStart:end);
            bkRightCC = sum(bkRight(currBkMapIx));

            p = realK/bkRightW;
            kmerRightPV = compute_PV(p, kmerCC, bkRightCC);            
            allPV(kmAllCount) = max(kmerLeftPV, kmerRightPV);
          end

          allStart(kmAllCount) = winStart(k,c);
          allStop(kmAllCount) = winStop(k,c);
          kmAllCount = kmAllCount + 1;
        end
      end
    end
  end

  [cleanPV, cleanStart, cleanStop] = ...
    remove_overlaps(allPV(1:kmAllCount-1), allStart(1:kmAllCount-1), allStop(1:kmAllCount-1));

  [junk, chrIx] = sort(cleanStart); 
  for i = 1:length(chrIx)
    fprintf(fpDataTemp, '%d\t%d\t%d\t%g\n', ...
            rowNumber, cleanStart(chrIx(i)), cleanStop(chrIx(i)), cleanPV(chrIx(i)));
  end


  % ----null model----
  nullRow = data(mappableIx);
  nullBkLen = localBkWidth;
  nullBkHalf = floor(nullBkLen / 2);
  nullFirstCol = 1 + flankSize;
  nullLastCol = validLen + flankSize;

  for s = 1:validLen-1
    rnd = floor(rand * (validLen-s))+1+s;
    tmp = nullRow(s);
    nullRow(s) = nullRow(rnd);
    nullRow(rnd)= tmp;
  end

  nullRow = [data(1:flankSize) ...
             nullRow ...
             data((end-flankSize+1):end)];
  nullMappable = [mappable(1:flankSize) ...
                  ones(1, validLen) ...
                  mappable((end-flankSize+1):end)];

  kmMaxNum = validLen * numK;
  allNullPV = repmat(-1, 1, kmMaxNum);
  allNullStart = repmat(-1, 1, kmMaxNum);
  allNullStop = repmat(-1, 1, kmMaxNum);
  kmNullCount = 1;
  
  for c = (1+flankSize):(validLen+flankSize)  
    nullBkStart = c - nullBkHalf;
    nullBkStop = c + nullBkHalf;

    nullBkMappable = nullMappable(nullBkStart:nullBkStop);
    
    if transformType == 0
      nullBk = nullRow(nullBkStart:nullBkStop);
    else
      nullBk = generate_input(nullRow(nullBkStart:nullBkStop), ...
                              nullBkMappable, ...
                              1, nullBkLen, transformType); 
    end

    if ~isOneSideBkg
      nullBkMapIx = find(nullBkMappable ~= 0);
      nullBkW = sum(nullBkMappable);
      nullBkCC = floor(sum(nullBk(nullBkMapIx)));
    end

    for k = minK:KStep:maxK

      if (winStart(k,c) >= nullFirstCol) && (winStop(k,c) <= nullLastCol)
        nullStart = winStart(k,c) - nullBkStart + 1;        
        nullStop = winStop(k,c) - nullBkStart + 1;

        kmNullCC = sum(nullBk(nullStart:nullStop));

        if ~isOneSideBkg
          p = k/nullBkW;
          kmNullPV = compute_PV(p, kmNullCC, nullBkCC);
          allNullPV(kmNullCount) = kmNullPV;
        else
          % left side background
          nullBkLeftW = sum(nullBkMappable(1:nullStop));
          nullBkMapIx = find(nullBkMappable(1:nullStop) ~= 0);
          nullBkLeft = nullBk(1:nullStop);
          nullBkLeftCC = sum(nullBkLeft(nullBkMapIx));

          p = k/nullBkLeftW;
          kmNullLeftPV = compute_PV(p, kmNullCC, nullBkLeftCC);

          % right side background
          nullBkRightW = sum(nullBkMappable(nullStart:end));
          nullBkMapIx = find(nullBkMappable(nullStart:end) ~= 0);
          nullBkRight = nullBk(nullStart:end);
          nullBkRightCC = sum(nullBkRight(nullBkMapIx));

          p = k/nullBkRightW;
          kmNullRightPV = compute_PV(p, kmNullCC, nullBkRightCC);          
          allNullPV(kmNullCount) = max(kmNullLeftPV, kmNullRightPV);
        end
       
        allNullStart(kmNullCount) = winStart(k,c);
        allNullStop(kmNullCount) = winStop(k,c);
        kmNullCount = kmNullCount + 1;
      end
    end
  end


  cleanNullPV = remove_overlaps(allNullPV(1:kmNullCount-1), allNullStart(1:kmNullCount-1), allNullStop(1:kmNullCount-1));

  for i = 1:length(cleanNullPV)
    fprintf(fpNullTemp, '%d\t%g\n', rowNumber, cleanNullPV(i));
  end

end