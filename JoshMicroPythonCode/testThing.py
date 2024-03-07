import time
# import thumby


book = open("C:/Users/Josh2/Desktop/book.txt", "r", encoding="utf8")
book = book.read()
book = book.replace('\n',' ')
# book = book.replace(' ','\n')

book = book.split(" ")
length = len(book)
lines = []
numChar = 72//6
i=0
while i+1<=length:
    word = book[i]
    if len(word) > numChar:
        lines = lines+[word+'\n']
        i+=1
    elif(i+1 == length):
        break
    else:
        line = word
        j = 1
        while(len(line+book[i+j])<=numChar):
            line = line+' '+book[i+j]
            if i+j+2 <= length:
                j += 1
            else:
                break              
        lines = lines+[line+'\n']
        i+=j

book = lines
line = 0
length = len(book)