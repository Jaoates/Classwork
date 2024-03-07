import time
import thumby

def abWait(s):
    thumby.display.fill(0)
    thumby.display.drawText("Press A/B",0,0,1)
    thumby.display.drawText(s,0,10,1)
    thumby.display.update()
    time.sleep(.1)
    while(not thumby.actionJustPressed()):
        pass


book = open("/Games/Reader/book.txt", "rt")
book = book.read()
book = book.replace('\n',' ')
# book = book.replace(' ','\n')

book = book.split(" ")
book = book+[" "," ","-File End-"," "]
length = len(book)
lines = []
numChar = thumby.display.width//6-1
i=0
while i+1<=length:
    word = book[i].replace(' ','')
    if len(word) > numChar:
        lines = lines+[word+'\n']
        i+=1
    elif(i+1 == length):
        break
    else:
        line = word
        j = 1
        while(len(line+book[i+j])<=numChar):
            line = line+' '+book[i+j].replace(' ','')
            if i+j+2 <= length:
                j += 1
            else:
                break              
        lines = lines+[line+'\n']
        i+=j

book = lines
line = 0
length = len(book)

abWait("to start")

while(1):
    thumby.display.fill(0)
    thumby.display.drawText(book[line],0,0,1)
    thumby.display.drawText(book[line+1],0,10,1)
    thumby.display.drawText(book[line+2],0,20,1)
    thumby.display.drawText(book[line+3],0,30,1)
    thumby.display.drawFilledRectangle(thumby.display.width-6*len(str(line))-1,thumby.display.height-8,6*len(str(line)),7,0)
    thumby.display.drawText(str(line),thumby.display.width-6*len(str(line)),thumby.display.height-7,1)
    thumby.display.update()

    jsize = 10
    if(thumby.buttonU.justPressed() and line>0):
        line -= 1
    elif(thumby.buttonD.justPressed() and line<length-4):
        line += 1
    elif(thumby.buttonR.justPressed() and line<length-4-jsize):
        line += jsize
    elif(thumby.buttonL.justPressed() and line-jsize>0):
        line -= jsize
    elif(thumby.actionJustPressed()):
        thumby.display.fill(0)
        thumby.display.drawText("Exiting...",0,0,1)
        thumby.display.update()
        time.sleep(1)
        break

