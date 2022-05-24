Program szukajacy motywu w sekwencji nukleotydowej. 

Kazda instancja sklada sie z dwoch plikow: 'inst.fasta' oraz 'inst.qual'. Pierwszy plik zawiera fragment 5 sekwencji nukleotydowych w formacie FASTA, drugi to ocena wiarygodnosci danego nukleotydu w danej sekwencji (kolejnosc danych fragmentow powinna pozostac taka sama w obu plikach).
Program usuwa nukleotydy ponizej danego progu wiarygodnosci a nastepnie szuka motywu (takich samcyh fragmnetow we wszystkich pieciu sekwencjach). Motyw szukany jest poprzez szukanie struktury kliko-podobnej w grafie. 

Program zrealizowany w ramach laboratorium "Algorytmy Kombinatoryczne w Bioinformatyce" prowadzonym przez prof. Martę Kasprzak. 