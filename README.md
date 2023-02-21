# Operační systémy 2. projekt synchronizace

# Popis Úlohy (Building H2O)
Molekuly vody jsou vytvářeny ze dvou atomů vodíku a jednoho atomu kyslíku. V systému jsou tři typy
procesů: (0) hlavní proces, (1) kyslík a (2) vodík. Po vytvoření procesů se procesy reprezentující
kyslíky a vodíky řadí do dvou front—jedna pro kyslíky a druhá pro vodíky. Ze začátku fronty vždy
vystoupí jeden kyslík a dva vodíky a vytvoří molekulu. V jednu chvíli je možné vytvářet pouze jednu
molekulu. Po jejím vytvoření je prostor uvolněn dalším atomům pro vytvoření další molekuly. Procesy,
které vytvořily molekulu následně končí. Ve chvíli, kdy již není k dispozici dostatek atomů kyslíku
nebo vodíku pro další molekulu (a ani žádné další již nebudou hlavním procesem vytvořeny) jsou
všechny zbývající atomy kyslíku a vodíku uvolněny z front a procesy jsou ukončeny.