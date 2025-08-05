class Find:
    def __init__(self, seq:str):
        self.seq = str(seq)
        
    def match(self, sequence:str):
        return (self.seq in str(sequence))
    
    def __and__(self, other):
        return And(self, other)
        
    def __or__(self, other):
        return Or(self, other)
    
    def __invert__(self):
        return Not(self)
    
    def __repr__(self):
        return f"Find \"{self.seq}\""
    
class And:
    def __init__(self, objA, objB):
        self.objA, self.objB = objA, objB
            
    def match(self, sequence:str):
        return self.objA.match(sequence) and self.objB.match(sequence)
    
    def __and__(self, other):
        return And(self, other)
    
    def __or__(self, other):
        return Or(self, other)
    
    def __invert__(self):
        return Not(self)
        
    def __repr__(self):
        return str(self.objA) + " & " + str(self.objB) + " "
        
class Or:
    def __init__(self, objA, objB):
        self.objA, self.objB = objA, objB
            
    def match(self, sequence:str):
        return self.objA.match(sequence) or self.objB.match(sequence)
    
    def __and__(self, other):
        return And(self, other)
    
    def __or__(self, other):
        return Or(self, other)
    
    def __invert__(self):
        return Not(self)
        
    def __repr__(self):
        return str(self.objA) + " | " + str(self.objB) + " "
    
class Not:
    def __init__(self, obj):
        self.obj = obj
    
    def match(self, sequence):
        return (not self.obj.match(sequence))
    
    def __and__(self, other):
        return And(self, other)
    
    def __or__(self, other):
        return Or(self, other)
    
    def __invert__(self):
        return Not(self)
    
    def __repr__(self):
        return "~" + str(self.obj)
    
if __name__ == "__main__":
    A = "ATCGG"
    B = "ACTGA"
    C = "ACGCG"
    seqs = ["ATCGGACTGAACGCG", "ATCGGAGTC", "ACTGAATAA"]
    query1 = Find(A) & Find(B) & Find(C)
    query2 = (Find(A) | Find(B)) & (~ Find(C))
    print([seq for seq in seqs if query1.match(seq)])
    print([seq for seq in seqs if query2.match(seq)])
    print(query1, query2)