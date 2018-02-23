
#Person 1
card_list =["1R","2R","3R","1G","2G","3G","1B","2B","3B"]

#Person 2
red_list = []

#Person 3
green_list = []

#Person 4
blue_list = []

#Person 9
number = 0

#Person 5
for card in card_list:
    
    #Person 6
    if card[-1] == "R":
        red_list = red_list + [card]
        
    #Person 7
    if card[-1] == "G":
        green_list = green_list + [card]
        
    #Person 8
    if card[-1] == "B":
        blue_list = blue_list + [card]
        
    #Person 10
    number = number + int(card[0])

#Person 11
print number
