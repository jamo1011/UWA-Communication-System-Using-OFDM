function [isiSymbol] = isi(symbol, previousSymbol, offset, amplitude)

isiSymbol=symbol;

for i=1:length(offset)
    
    interference = [previousSymbol(length(previousSymbol)-offset(i):length(previousSymbol)),zeros(1,length(previousSymbol)-offset(i)-1)].*amplitude(i);
    isiSymbol = isiSymbol+interference;
    
end


end