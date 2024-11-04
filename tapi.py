import pandas as pd
import numpy as np
import math

# קריאת קובץ אקסל לדאטה פריים
file_path = input("Enter the file path: ")
df = pd.read_excel(file_path, header=None, names=['charges'])
df = df.iloc[1:, :]  #חילוץ העמודה הראשונה כולל כותרת

def path_excel():
    path=file_path
    return path

def c(df): # פונקציה שמחזירה את מקדמי העונתיות
    new_df = pd.DataFrame(columns=range(1, 13),index=range(1, 7))  # # דאטה פריים בתוספת שורה לחישובי ממוצע ומקדמי עונתיות
    tempdf = pd.DataFrame(columns=range(1, 13), index=range(1, 6))  # דאטה פריים מקורי שעליה נעשה נטרול עונתיות
    # מילוי דאטה פריימים נקיים (מתקבלים כערכים בפונקציה) מנתוני האקסל (df)
    for i in range(1, 6):
        for j in range(1, 13):
            new_df.loc[i, j] = df.iloc[(i-1)*12 + j-1, 0]
            tempdf.loc[i,j]=df.iloc[(i-1)*12 + j-1, 0]
    new_df['average']=new_df.mean(axis=1) # הוספת טור עם ממוצע כל מחזור

    for i in range(1, 6):
        for j in range(1, 13):
            new_df.loc[i,j]=new_df.loc[i,j]/new_df.loc[i,'average']#חלוקה כל ערך בממוצע המחזור שלו-ניטרול עונתיות
    for i in range(1,13):
        new_df.loc[6,i]=new_df[i].mean()#מוסיף בשורה האחרונה את ממוצע כל עונה- מקדמי העונתיות3
    x=[] # רשימה שבה ישמרו מקדמי העונתיות
    c=new_df.loc[6,1:12]#מערך של מקדמי העונתיות
    for i in range(1, 13):
        x.append(float(c.loc[i]))
    return x

def if_lyniar_fixed_exponntial(df):
    a=[]
    x=[]
    filtered=[]
    for i in range(1,31):
        a.append(abs(float(df.loc[i])-float(df.loc[i+1])))
    for i in range(len(a)-1):
        x.append(abs(a[i+1]-a[i]))
    UCL = np.percentile(x, 90)
    LCL = np.percentile(x, 10)
    for i in x:
        if i <=UCL and i >=LCL:
            filtered.append(i)
    avg=np.mean(filtered)
    if avg<1:
        return True
    return False

def decomposition(c): # פונקציה שמבצעת דה- קומפוזיציה
    c1=pd.DataFrame({'1': c}) # מכניסים את רשימת מקדמי העונתיות למערך
    c1.index += 1
    tempdf = pd.DataFrame(columns=range(1, 13), index=range(1, 6))  # דאטה פריים מקורי שעליה נעשה נטרול עונתיות
    for i in range(1, 6): # 2 לולאות שמכניסות לtempdf את הנתונים
        for j in range(1, 13):
            tempdf.loc[i,j]=df.iloc[(i-1)*12 + j-1, 0]
    season_average = tempdf.div(c1["1"], axis=1)  # חלוקת כל עונה במקדם העונתיות המתאים
    season_average['average'] = season_average.mean(axis=1)  # חישוב ממוצע של מחזור לאחר נטרול עונתיות
    return pd.DataFrame({'1': season_average["average"]})

def reg2(season_average): #פונקציה שמקבלת מערך ומחזירה את אומדי הרגרסיה a b במודל עונתי ריבועי
    n=len(season_average)
    Dt=season_average.sum()
    t2=0
    t4=0
    Dt2=0
    for i in range(1,n+1):
        t2+=i**2
        t4+=i**4
        Dt2+=season_average.loc[i]*(i**2)
    b=(n*Dt2-Dt*t2)/(n*t4-t2**2)
    a=(Dt-b*t2)/n
    return (a,b)

def season_square_mse(season_average): # חישוב MSE למודל עונתי ריבועי
    mse=square_mse(season_average)
    return mse

def f2(season_average,c): # פונקציה שמייצרת תחזית ל12 חודשים הבאים במודל עונתי ריבועי
    x=[]
    ab = reg2(season_average)
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range(len(c)):
        x.append(float(c[i])*(ab[0]+ab[1]*6**2))
    for i in range(61,73):
        ft.loc[i]=float(x[i-61])
    return (ft,"FT = c_t * ("+str(float(ab[0]))+str(float(ab[1]))+" * 6^2)")

def reg1(season_average):#פונקציה שמקבלת מערך ומחזירה את אומדי הרגרסיה a b במודל עונתי ליניארי
    n=len(season_average)
    Dt=season_average.sum()
    t=0
    t2=0
    tDt=0
    for i in range(1,n+1):
        t+=i
        t2+=i**2
        tDt+=season_average.loc[i]*i
    b=(n*tDt-Dt*t)/(n*t2-t**2)
    a=(Dt-b*t)/n
    return (a,b)

def season_lyniar_reg_mse(season_average):
    mse=lyniar_reg_mse(season_average)
    return mse

def season_lyniar_holt_alpha_beta_mse(season_average):
    mse=best_alpha_beta_mse(season_average)
    return mse

def season_lyniar_brown_alpha_mse(season_average):
    mse=brown_best_alpha_mse(season_average)
    return mse

def season_lyniar_double_moving_k_mse(season_average):
    mse=double_best_k_mse(season_average)
    return mse

def f1(season_average,c):#תחזית במידה והמודל עונתי ליניארי
    ab=reg1(season_average)
    x=[]
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range(len(c)):
        x.append(c[i]*(ab[0]+ab[1]*6))
    for i in range (61,73):
        ft.loc[i] = float(x[i - 61])
    return (ft,"FT = c_t * ("+str(float(ab[0]))+"+ ("+str(float(ab[1]))+") * 6)")

def season_avg_fixed_mse(season_average):
    mse=fixed_avg_mse(season_average)
    return mse

def season_moving_avg_k_mse(season_average):
    mse=best_k_mse(season_average)
    return mse

def season_exponntial_alpha_mse(season_average):
    mse=best_alpha_mse(season_average)
    return mse

def reg(season_average): #פונקציה שמקבלת מערך ומחזירה את אומד a במודל עונתי קבוע בזמן
    n=len(season_average)
    Dt=season_average.sum()
    a=Dt/n
    return a

def f(season_average,c):# תחזית במידה והמודל עונתי קבוע בזמן
    a=reg(season_average)
    x=[]
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range(len(c)):
        x.append(c[i]*a)
    for i in range (61,73):
        ft.loc[i] = float(x[i - 61])
    return (ft,"FT = c_t * "+str(float(a)))

#----------------------------------מודל קבוע בזמן------------------------------

def fixed_avg(df):  # חיזוי לכל התקופות במודל קבוע בזמן עם ממוצע פשוט
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    f = df.mean()
    for i in range (61,73): #לולאה שמייצרת חיזוי ל12 תקופות הבאות
        ft.loc[i]=float(f)
    return (ft ,"\nFt = "+str(f))

def fixed_avg_mse(df): # מודל קבוע בזמן ממוצע פשוט חישוב MSE
    ft =[]
    sum=0
    for i in range(len(df)-1): #חישוב חיזוי בדיעבד
        sum+=df.loc[i+1]
        avg=sum/(i+1)
        ft.append(float(avg))
    e =[] #רשימה שבה יהיו הטעויות בריבוע
    for i in range(len(df)-1): # לולאה שמחשבת את הטעות ומכניסה לרשימה
        e.append((ft[i]-float(df.loc[i+2]))**2)
    sum=0
    for i in e:
        sum+=i
    mse=sum/(len(df)-1)
    return mse

def best_k_mse(df): # פונקציה שבודקת מה הK שעבורו מתקבל הMSE הנמוך ביותר לממוצע נע , מחזירה K וMSE
    mini_mse = 100000000 #אתחול הMSE המינימלי
    best_k = 60

    for k in range(1, len(df)): # לולאה שעוברת על כל הK האפשריים
        a = []
        for i in range(1, len(df) - (k - 1)): # לולאה שמבצעת חיזוי בדיעבר על כל תקופה לפי ממוצע נע עם הK הרלוונטי
            sum = 0
            for j in range(i, i + k):
                sum += df.loc[j]
            avg = sum / k
            a.append(avg)

        e = []  # רשימה שבה יהיו הטעויות בריבוע
        for i in range(len(df) - k): # לולאה שמחשבת את הטעויות
            e.append((a[i] - float(df.loc[i + k + 1])) ** 2)
        sum1=0
        for i in e:
            sum1+=i
        mse = sum1 / (len(df) - k)

        if float(mse) < mini_mse: # לוקח את הMSE המינימלי ושמירת הK שעבורו זה מתקיים
            mini_mse = float(mse)
            best_k = k
    return (best_k,mini_mse)

def fixed_moving_k(df,k):  # חיזוי לכל התקופה במודל קבוע בזמן בעזרת ממוצע נע עם K מהפונקציה שמחזירה K הכי טוב
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    sum = 0
    if k !=60: #תנאי שמוודא שזה לא מקרה של ממוצע פשוט
        for i in range(len(df) - k + 1, len(df) + 1): #לולאה שמבצעת חיזוי עם K נתונים אחרונים
            sum += float(df.loc[i])
        f=sum / k
        for i in range (61,73):
            ft.loc[i]=f
        return (ft,"\n Ft="+str(f))
    return fixed_avg(df)

def best_alpha_mse(df): # פונקציה שבודקת מה האלפא שעבורה מתקבל הMSE הנמוך ביותר להחלקה אקספוננציאלית,מחזירה אלפא וMSE
    best_alpha = 1
    mini_mse=1000000000
    for alpha in range(10,31,1): # לולאה שעוברת על כל האלפות האפשריות
        e = [] # אתחול רשימה שבה יהי הטעויות בריבוע
        alpha/=100
        a = [] #אתחול רשימה שבה יהיה החיזוי
        a.append(float(df.loc[1]))
        for i in range(2, len(df)):  # חישוב חיזוי בדיעבד לפי החלקה אקספוננציאלית לכל אלפא
            a.append(alpha * float(df.loc[i]) + (1-alpha) * a[i - 2])
        for i in range(2, len(a) + 1): # מחשב את הטעות ומוסיף אותה לרשימת הטעויות
            e.append((a[i - 2] - float(df.loc[i])) ** 2)
        mse = (sum(e)) / (len(df)-1)
        if float(mse) < mini_mse: # לוקח את הMSE המינימלי ושמירת האלפא שעבורה זה מתקיים
            mini_mse = float(mse)
            best_alpha = alpha
    return (best_alpha,mini_mse)

def fixed_exponntial_smooting(df,alpha): # פונקציה שמחזירה חיזוי לפי החלקה אקספונינציאלית עם אלפא שיצאה בפונקציה של האלפא הכי טובה
    ft = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    a = [] # רשימה שבה יהיה החיזוי
    a.append(float(df.loc[1]))
    for i in range(1, len(df)):  # לולאה שמחשבת עבור כל תקופה את החיזוי שלה לפי החלקה אקספוננצילית
        a.append(alpha * float(df.iloc[i]) + (1-alpha) * a[i - 1])
    for i in range(61,73):
        ft.loc[i]=a[len(a) - 1]
    return (ft,"\nFt = "+str(a[len(a)-1]))

#----------------------------------מודל ליניארי------------------------------

def lyniar_reg_mse(df): # פונקציה שמחשבת MSE לחיזוי עפ ריגרסיה ליניארית
    ft = []
    for i in range(2, len(df)):  # לולאה שמחשבת חיזוי בדיעבד
        f=df.iloc[:i,:] # חיתוך DT כדי לחשב חיזוי לתקופה הבאה
        ab=reg1(f)
        temp=float(ab[0])+float(ab[1])*(i+1)
        ft.append(temp)
    e = []  # מערך שבו יהיו הטעויות
    for i in range(3, len(df) + 1): # לולאה שמחשבת את הטעויות בריבוע
        e.append ((float(ft[i-3]) - float(df.loc[i])) ** 2)
    mse = (sum(e)) / (len(e)) # חישוב MSE
    return float(mse)

def lyniar_reg (df): # פונקציה שמחזירה חיזוי למודל ליניארי לפי רגרסיה ליניארית
    ab=reg1(df) # זימון הפונקציה שמחזירה את האומדים לA,B
    f = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range (61,73):
        f.loc[i]=float(ab[0])+float(ab[1])*i
    return (f,"Ft = "+str(float(ab[0]))+" + ("+str(float(ab[1]))+" * t)")

def double_best_k_mse(df): # פונקציה שמחשבת מהו הK שעבורו מתקבל הMSE הנמוך ביותר בחישוב ממוצע נע כפול, מחזירה K וMSE

    chek_df = pd.read_excel(path_excel(), header=None, names=['charges'])
    chek_df = chek_df.iloc[1:, :]  # קריאה מהקובץ המקורי להשוואת הדאטה שקיבלנו בפונקציה על מנת לדעת אם הנתונים מותמרים לוגריתמית
    flag = False#דגל איתות לפונקציה שאומר אם הנתונים שקיבלנו מותמרים או לא
    if len(df) == 60:
        for l in range(1, len(df) + 1):
            if float(chek_df.loc[l]) != float(df.loc[l]):#אם לפחות אחד מהנתונים שקיבלנו שונה מנתוני האקסל המקוריים האיתות יופעל (סימן שהנתונים מותמרים)
                flag = True
    mini_mse = 1000000000
    best_k = 30 # מקסימום K יכול להיות 2K-1=60
    e=[] #רשימה שבה יכנסו הטעויות בריבוע
    for k in range(2,31):  # לולאה שרצה על כל הK האפשריים
        mt1= []
        mt2=[]
        e = []  # רשימה שבה יכנסו הטעויות בריבוע
        a = [] # רשימה שבה יהיה החיזוי בדיעבד
        for i in range(1, len(df) - (k - 2)):  # לולאה שרצה כל DT ובונה את MT1
            temp=df.loc[i:i+k-1,:]
            mt1.append(float(temp.mean()))
        for j in range(len(mt1)-k+1): #לולאה שרצה על MT1 ובונה את MT2
            temp2=mt1[j:k+j]
            mt2.append(sum(temp2)/len(temp2))
        for p in range(len(mt2)-1): #לולאה שרצה באורך של MT2 (פחות האיבר האחרון) ובונה חיזוי בדיעבד לכל תקופה רלוונטית
            dt=2*mt1[p+k-1]-mt2[p]
            bt=(2/(k-1))*(mt1[p+k-1]-mt2[p])
            a.append(dt+bt)
        if flag==True:
            y=pd.DataFrame(columns=range(1, 2), index=range(1, 61))
            for z in range (1,len(df)+1):
                y.loc[z]=math.e**float(df.loc[z])
            for s in range (len(a)):
                a[s] = math.e**a[s]#המרת החיזוי בדיעבד לרמה תואמת של נתוני האמת שקיבלו
            for l in range(len(a)):  # לולאה שמוסיפה לרשימה e את השגיאות בריבוע
                e.append(float((a[l] - y.loc[2 * k + l]) ** 2))
            if len(e) != 0:
                mse = sum(e) / len(e)  # חישוב MSE עבור כל K
                if mse < mini_mse:  # לוקח את הMSE המינימלי ושמירת הK שעבורו זה מתקיים
                    mini_mse = mse
                    best_k = k
            return (best_k, mini_mse)
        for l in range (len(a)): #לולאה שמוסיפה לרשימה e את השגיאות בריבוע
            e.append(float((a[l]-df.loc[2*k+l])**2))
        if len(e)!=0:
            mse=sum(e)/len(e)# חישוב MSE עבור כל K
            if mse<mini_mse: #לוקח את הMSE המינימלי ושמירת הK שעבורו זה מתקיים
                mini_mse=mse
                best_k=k
    return (best_k,mini_mse)

def lyniar_double_moving_k(df,k): #חיזוי במודל ליניארי בעזרת ממוצע נע כפול
    mt1 = []
    mt2 = []
    f= pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range(1, len(df) - (k - 2)):  # לולאה שרצה כל DT ובונה את MT1
        temp = df.loc[i:i + k - 1, :]
        mt1.append(float(temp.mean()))
    for j in range(len(mt1) - k + 1):  # לולאה שרצה על MT1 ובונה את MT2
        temp2 = mt1[j:k + j]
        mt2.append(sum(temp2) / len(temp2))
    mt1_60=mt1[len(mt1)-1] #האיבר האחרון בMT1
    mt2_60=mt2[len(mt2)-1] # האיבר האחרון בMT2
    D60=2*mt1_60-mt2_60 #אומדן לD האחרון
    b60=mt1_60-mt2_60 # אומדן לשיפוע הקו
    for i in range(1,13): # יצירת החיזוי בעזרת מערכת צירים דינאמית
        f.loc[i+60]=D60+b60*i
    return (f,"Ft_60+τ = "+str(D60)+" + ("+str(b60)+" * τ)",D60,b60)

def best_alpha_beta_mse(df): # פונקציה שמחשבת את האלפא ובטא שעבורם מתקבל הMSE הנמוך ביותר בחישוב החלקה אקספוננציאלית (הולט), מחזיר אלפא בטא וMSE
    mini_mse=10000000
    best_alpha=0.3
    best_beta=0.3
    chek_df = pd.read_excel(path_excel(), header=None, names=['charges'])
    chek_df = chek_df.iloc[1:, :]  # קריאה מהקובץ המקורי להשוואת הדאטה שקיבלנו בפונקציה על מנת לדעת אם הנתונים מותמרים לוגריתמית
    flag = False #דגל איתות לפונקציה שאומר אם הנתונים שקיבלנו מותמרים או לא
    if len(df) == 60:
        for l in range(1, len(df) + 1):
            if float(chek_df.loc[l]) != float(df.loc[l]):
                flag = True
    for alpha in range(1,4): # לולאה שעוברת על כל האלפא האפשריות
        alpha /= 10
        for beta in range(1,4): # לולאה שעוברת על כל הבטא האפשריות
            e = [] # רשימה שבה יהיו כל הטעויות
            a = [] # רשימה שבה יהיה החיזוי
            beta/= 10
            for i in range(2, len(df)):
                temp = df.iloc[:i, :] #חיתוך הDT לביצוע רגרסיה לפי תקופות
                ab = reg1(temp) # זימון הפונקציה שמוציאה אומדים לA,B לפי רגרסיה
                d0 = ab[0] # אתחול D0
                b0 = ab[1] # אתחול B0
                for t in range (len(temp)): # לולאה שמבצעת חיזוי בדיעבד לתקופה הבאה
                    di=float(alpha*float(df.loc[t+1])+(1-alpha)*(d0+b0))
                    bi=float(beta*(di-d0)+(1-beta)*b0)
                    d0=di
                    b0=bi
                a.append(di+bi) # נוסחת החיזוי לתקופה אחת קדימה (t=1)
            if flag == True:
                y=pd.DataFrame(columns=range(1, 2), index=range(1, 61))
                for z in range (1,len(df)+1):
                    y.loc[z]=math.e**float(df.loc[z])
                for s in range(len(a)):
                    a[s] =math.e**a[s]#המרת החיזוי בדיעבד לרמה תואמת של נתוני האמת שקיבלו
                for i in range(len(y) - 2):  # לולאה שמחשבת את הטעות בריבוע ומוסיפה אותה לרשימה
                    e.append(float((a[i] - y.loc[i + 3]) ** 2))
                mse = (sum(e) / len(e))  # חישוב MSE
                if mse < mini_mse:  # לוקח את הMSE המינימלי ושמירת האלפא ובטא שעבורן זה מתקיים
                    mini_mse = mse
                    best_alpha = alpha
                    best_beta = beta
                return (best_alpha, best_beta, mini_mse)
            for i in range(len(df)-2): # לולאה שמחשבת את הטעות בריבוע ומוסיפה אותה לרשימה
                e.append(float((a[i]-df.loc[i+3])**2))
            mse=(sum(e)/len(e)) # חישוב MSE
            if mse<mini_mse: #לוקח את הMSE המינימלי ושמירת האלפא ובטא שעבורן זה מתקיים
                mini_mse=mse
                best_alpha=alpha
                best_beta=beta
    return (best_alpha,best_beta,mini_mse)

def lyniar_exponntial_smooting(df,alpha,beta):# פונקציה שמבצעת חיזוי לפי החלקה אקספוננציאלית עם אלפא ובטא (מודל הולט)
    ab=reg1(df) # זימון הפונקציה שמחזירה את האומדים לA,B לפי רגרסיה
    d0=ab[0] # אתחול D0
    b0=ab[1] # אתחול B0
    f = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range(1, len(df)+1): # לולאה שמייצרת אומדים לD וB לכל תקופה
        dt=alpha*float(df.loc[i])+(1-alpha)*(d0+b0)
        bt=beta*(dt-d0)+(1-beta)*b0
        d0=dt
        b0=bt
    for i in range(1,13): # יצירת החיזוי בעזרת מערכת צירים דינאמית
        f.loc[60+i]=float(d0)+float(b0)*i
    return (f,"Ft_60+τ = "+str(float(d0))+" + ("+str(float(b0))+" * τ)",d0,b0)

def brown_best_alpha_mse(df): # פונקציה שמחשבת את האלפא שעבורה מתקבל הMSE הנמוך ביותר בחישוב החלקה אקספוננציאלית כפולה (בראון), מחזיר אלפא וMSE
    best_alpha=0.3
    mini_mse=100000000000000
    chek_df = pd.read_excel(path_excel(), header=None, names=['charges'])
    chek_df = chek_df.iloc[1:, :]  # קריאה מהקובץ המקורי להשוואת הדאטה שקיבלנו בפונקציה על מנת לדעת אם הנתונים מותמרים לוגריתמית
    flag = False#דגל איתות לפונקציה שאומר אם הנתונים שקיבלנו מותמרים או לא
    if len(df) == 60:
        for l in range(1, len(df) + 1):
            if float(chek_df.loc[l]) != float(df.loc[l]):
                flag = True
    for alpha in range(10,31,1): # לולאה שעוברת על כל האלפא האפשריות
        alpha/=100
        a=[] # רשימה שבה יהיה החיזוי
        e=[] # רשימה שבה יהיו הטעויות
        for i in range(2,len(df)): # לולאה שמחשבת S01 וS02 וחיזוי בדיעבד
            temp=df.iloc[:i,:] # חיתוך DF לפי התקופות
            ab=reg1(temp) # זימון הפונקציה שמייצרת אומד לA,B ברגרסיה
            s01=ab[0]-((1-alpha)/alpha)*ab[1] # אתחול S01
            s02=ab[0]-2*((1-alpha)/alpha)*ab[1] # אתחול S02
            for t in range(1,len(temp)+1): # לוללאה שמחשבת ST1, ST2
                st1=alpha*float(df.loc[t])+(1-alpha)*s01
                st2=alpha*st1+(1-alpha)*s02
                s01=st1
                s02=st2
            dt=2*s01-s02 # חישוב DT לפי ST1,ST2 האחרונים
            bt=(alpha/(1-alpha))*(s01-s02) # חישוב BT לפי ST1,ST2 האחרונים
            a.append(float(dt+bt)) # חיזוי בדיעבד לתקופה הבאה
        if flag==True:
            y=pd.DataFrame(columns=range(1, 2), index=range(1, 61))
            for z in range (1,len(df)+1):
                y.loc[z]=math.e**float(df.loc[z])
            for s in range (len(a)):
                a[s] = math.e**a[s]#המרת החיזוי בדיעבד לרמה תואמת של נתוני האמת שקיבלו
            for j in range(len(y) - 2):  # לולאה שמחשבת טעויות בריבוע
                e.append(float((a[j] - y.loc[j + 3]) ** 2))
            mse = (sum(e) / len(e))  # חישוב MSE
            if mse < mini_mse:  # לוקח את הMSE המינימלי ושמירת האלפא שעבורה זה מתקיים
                mini_mse = mse
                best_alpha = alpha
            return (best_alpha, mini_mse)
        for j in range(len(df) - 2): # לולאה שמחשבת טעויות בריבוע
            e.append(float((a[j] - df.loc[j + 3]) ** 2))
        mse = (sum(e) / len(e)) # חישוב MSE
        if mse < mini_mse: #לוקח את הMSE המינימלי ושמירת האלפא שעבורה זה מתקיים
            mini_mse = mse
            best_alpha = alpha
    return (best_alpha,mini_mse)

def lyniar_double_exponntial_smooting(df,alpha):#פונקציה שמחשבת חיזוי עבור החלקה אקספוננציאלית כפולה לפי האלפא מהפונקציה של האלפא הכי טובה
    ab=reg1(df) #זימון פונקציה שמייצרת אומדים לA,B לפי רגרסיה
    a=ab[0] # אתחול A
    b=ab[1] # אתחול B
    s01=a-((1-alpha)/alpha)*b # אתחול S01
    s02=a-2*((1-alpha)/alpha)*b # אתחול S02
    f = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    st1= []
    st2 = []
    for i in range(1, len(df) + 1): # פונקציה שמייצרת מערך של ST1
        st1.append(alpha*float(df.loc[i])+(1-alpha)*float(s01))
        s01=st1[i-1]
    for i in range(1, len(df) + 1): # פונקציה שמייצרת מערך של ST2
        st2.append(alpha*st1[i-1]+(1-alpha)*float(s02))
        s02=st2[i-1]

    d60=2*float(st1[len(st1)-1])-float(st2[len(st2)-1]) # חישוב הD האחרון לטובת החיזוי
    b60=(alpha/(1-alpha))*(float(st1[len(st1)-1])-float(st2[len(st2)-1])) # חישוב הB האחרון לטובת החיזוי

    for i in range(1,13): # יצירת החיזוי בעזרת מערכת צירים דינאמית
        f.loc[i+60]=float(d60)+(b60)*i
    return (f,"Ft_60+τ = "+str(d60)+" + ("+str(b60)+"*τ)",d60,b60)

#----------------------------------מודל ריבועי------------------------------

def square_f(df): #פונקציה שמייצרת חיזוי עבור מודל ריבועי בעזרת רגרסיה
    ab=reg2(df) # זימון פונקציה שמוציאה אומדים לA,B
    f = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    for i in range (61,73):
        f.loc[i]=float(ab[0])+float(ab[1])*(i**2)
    return (f,"Ft = "+str(float(ab[0]))+" + ("+str(float(ab[1]))+" * t^2)")

def square_mse(df): # פונקציה שמוצאת MSE לרגרסיה במודל ריבועי
    ft = []
    for i in range(2, len(df)):  # לולאה שמחשבת חיזוי בדיעבד
        f = df.iloc[:i, :] # חיתוך DT
        ab = reg2(f) # זימון פונקציה שמוציאה אומדים לA,B
        temp = float(ab[0]) + float(ab[1]) * (i + 1)**2
        ft.append(temp)
    e = []  # רשימה שבה יהיה הטעויות בריבוע
    for i in range(3, len(df) + 1): # לולאה שמחשבת את הטעות
        e.append((float(ft[i-3]) - float(df.loc[i])) ** 2)
    mse = np.mean(e)
    return float(mse)

#----------------------------------מודל מעריכי------------------------------

def e_reg(df): # פונקציה שמייצרת חיזוי בעזרת התמרה ורגרסיה ליניארית
    f = pd.DataFrame(columns=range(1, 2), index=range(61, 73))
    ft = pd.DataFrame(columns=range(1, 2), index=range(1, 61)) # מערך שבו יהיו הערכים המותמרים של נתוני החודשים
    for i in range (1,len(df)+1): # לולאה שעושה התמרה על DT
        ft.loc[i]=math.log(df.loc[i],math.e)
    ab=reg1(ft) # זימון של פונקציה שעושה רגרסיה ומחזירה אומדים לA,B
    b=ab[1]
    a=math.e**ab[0] # התמרה
    for i in range(61,73): # לולאה שמייצרת חיזוי ל12 תקופות הבאות
        f.loc[i]=float(a)*math.e**(float(b)*i)
    return (f,"FT = "+str(float(a))+" * e^("+str(float(b))+"*t)")

def e_reg_mse(df): # פונקציה שמוצאת MSE לחזוי בעזר רגרסיה על מודל מעריכי
    ft=[]
    e=[] # רשימה שבה יהיו הטעויות בריבוע
    for i in range(2,len(df)):
        x=[]
        temp=df.iloc[:i,:] # חיתוך של DT
        for j in range(1,len(temp)+1): #לולאה שמכניסה את DT החתוך לתוך רשימה
            x.append(float(temp.loc[j]))
        for k in range(len(x)): # לולאה שמבצעת התמרה על DT החתוך
            x[k]=math.log(x[k],math.e)
        x=pd.DataFrame(x,index=range(1, len(x) + 1))
        ab=reg1(x) # זימון פונקציה שמוציאה אומדים לA,B
        b=ab[1]
        a=math.e**ab[0]
        ft.append(float(a)*math.e**(float(b)*(i+1))) # חיזוי לתקופה הבאה
    for j in range(len(ft)): # חישוב הטעות בריבוע והכנסה לרשימה
        e.append((float(ft[j])-float(df.loc[j+3]))**2)
    mse=sum(e)/len(e) # חישוב MSE
    return mse

def e_double_moving_f_ft_k_mse(df): # פונקציה שמחזירה חיזוי לפי ממוצע נע כפול על ה K הכי טוב
    ft = pd.DataFrame(columns=range(1, 2), index=range(1, 61)) # מערך שבו יהיו הערכים המותמרים של נתוני החודשים
    for i in range (1,len(df)+1): # לולאה שמבצעת התמרה על הDT
        ft.loc[i]=math.log(df.loc[i],math.e)
    k=double_best_k_mse(ft) #  זימון הפונקציה שמוצאת את הK הכי טוב וה MSE
    f=lyniar_double_moving_k(ft,k[0]) # זימון פונקציה שמייצרת חיזוי על הDT עם ההתמרה
    new_f=math.e**f[0] # חזרה למודל מעריכי
    return (new_f,"Ft_60+τ = "+str(math.e**f[2])+" * e^("+str(f[3])+"*τ)",k[0],k[1])

def e_holt_f_ft_alpha_beta_mse(df):# פונקציה שמחזירה חיזוי לפי מודל הולט עם האלפא ובטא הכי טובות
    ft = pd.DataFrame(columns=range(1, 2), index=range(1, 61)) # מערך שבו יהיו הערכים המותמרים של נתוני החודשים
    for i in range (1,len(df)+1):  # לולאה שמבצעת התמרה על הDT
        ft.loc[i]=math.log(df.loc[i],math.e)
    alpha_beta=best_alpha_beta_mse(ft)#  זימון הפונקציה שמוצאת את האלפא ובטא הכי טובות וה MSE
    f=lyniar_exponntial_smooting(ft,alpha_beta[0],alpha_beta[1]) # זימון פונקציה שמייצרת חיזוי על הDT עם ההתמרה
    new_f=math.e**f[0] # חזרה למודל מעריכי
    return (new_f,"Ft_60+τ = "+str(math.e**float(f[2]))+" * e^("+str(float(f[3]))+"*τ)",alpha_beta[0],alpha_beta[1],alpha_beta[2])

def e_brown_f_ft_alpha_mse(df): # פונקציה שמחזירה חיזוי לפי מודל בראון עם האלפא  הכי טובה
    ft = pd.DataFrame(columns=range(1, 2), index=range(1, 61)) # מערך שבו יהיו הערכים המותמרים של נתוני החודשים
    for i in range (1,len(df)+1): # לולאה שמבצעת התמרה על הDT
        ft.loc[i]=math.log(df.loc[i],math.e)
    alpha=brown_best_alpha_mse(ft) #  זימון הפונקציה שמוצאת את האלפא הכי טובה וה MSE
    f=lyniar_double_exponntial_smooting(ft,alpha[0]) # זימון פונקציה שמייצרת חיזוי על הDT עם ההתמרה
    new_f=math.e**f[0] # חזרה למודל מעריכי
    return (new_f,"Ft_60+τ = "+str(math.e**float(f[2]))+" * e^("+str(f[3])+"*τ)",alpha[0],alpha[1])

#----------------------------------אופטימליים ------------------------------

def o(df):
    print("calculate...")
    ab = reg1(df)
    mse = {}
    if float(ab[1]) == 0:  # סינון מודלים קבועים בזמן באופן מובהק (ללא סטייה)
        k = best_k_mse(df)
        alpha = best_alpha_mse(df)
        regular_avg = fixed_avg_mse(df)
        mse[regular_avg] = "regular avg"
        mse[k[1]] = "mooving avg"
        mse[alpha[1]] = "exponntial smooting"
        best_mse = min(mse)
        if (mse[best_mse] == "regular avg") or (mse[best_mse] == "mooving avg" and k[1] == regular_avg):
            regular_f = fixed_avg(df)
            print("The model is fixed\n",regular_f[1],"\nThe forecast is made using regular avg","\nThe forecast is:\n",regular_f[0], "\nThe mse is:", best_mse)
        elif mse[best_mse] == "mooving avg":
            moving_f = fixed_moving_k(df, k[0])
            print("The model is fixed\n",moving_f[1],"\nThe forecast is made using ", mse[best_mse], "with k=", k[0],"\nThe forecas is:\n", moving_f[0], "\nThe mse is:", best_mse)
        else:
            smooting_f = fixed_exponntial_smooting(df, alpha[0])
            print("The model is fixed\n",smooting_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", alpha[0],"\nThe forecas is:\n", smooting_f[0], "\nThe mse is:", best_mse)
    elif if_lyniar_fixed_exponntial(df)==True:
        k = best_k_mse(df)
        alpha = best_alpha_mse(df)
        regular_avg = fixed_avg_mse(df)
        reg_mse = lyniar_reg_mse(df)
        holt = best_alpha_beta_mse(df)
        print("A few more seconds and you have the forecast!!\n")
        brown = brown_best_alpha_mse(df)
        double_moving = double_best_k_mse(df)
        reg_e = e_reg_mse(df)
        e_double = e_double_moving_f_ft_k_mse(df)
        e_brown = e_brown_f_ft_alpha_mse(df)
        e_holt = e_holt_f_ft_alpha_beta_mse(df)
        mse[reg_e] = "exponntial- regresia"
        mse[e_double[3]] = "exponntial- double moving"
        mse[e_brown[3]] = "exponntial- model brown"
        mse[e_holt[4]] = "exponntial- model holt"
        mse[reg_mse] = "rgresia"
        mse[holt[2]] = "model holt"
        mse[brown[1]] = "model brown"
        mse[double_moving[1]] = "double moving avg"
        mse[regular_avg] = "regular avg"
        mse[k[1]] = "mooving avg"
        mse[alpha[1]] = "exponntial smooting"
        best_mse=min(mse)
        if (mse[best_mse] == "regular avg") or (mse[best_mse] == "mooving avg" and k[1] == regular_avg):
            regular_f = fixed_avg(df)
            print("The model is fixed\n",regular_f[1],"\nThe forecast is made using regular avg", "\nThe forecast is:\n", regular_f[0]
                  , "\nThe mse is:", best_mse)
        elif mse[best_mse] == "mooving avg":
            moving_f = fixed_moving_k(df, k[0])
            print("The model is fixed\n",moving_f[1],"\nThe forecast is made using ", mse[best_mse], "with k=", k[0],
                  "\nThe forecast is:\n", moving_f[0], "\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial smooting":
            smooting_f = fixed_exponntial_smooting(df, alpha[0])
            print("The model is fixed\n",smooting_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", alpha[0],
                  "\nThe forecast is:\n", smooting_f[0], "\nThe mse is:", best_mse)
        elif mse[best_mse] == "rgresia":
            reg1_f = lyniar_reg(df)
            print("The model is lyniar\n",reg1_f[1],"\nThe forecast is made using ", mse[best_mse], "\nThe forecast is:\n", reg1_f[0]
                  , "\nThe mse is:", best_mse)
        elif mse[best_mse] == "model holt":
            holt_f = lyniar_exponntial_smooting(df, holt[0], holt[1])
            print("The model is lyniar\n",holt_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", holt[0],
                  "\nand beta=", holt[1],"\nThe forecast is:\n",
                  holt_f[0], "\nThe mse is:", best_mse)
        elif mse[best_mse] == "model brown":
            brown_f = lyniar_double_exponntial_smooting(df, brown[0])
            print("The model is lyniar\n",brown_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", brown[0],
                  "\nThe forecast is:\n",brown_f[0], "\n The mse is:", best_mse)
        elif mse[best_mse] == "double moving avg":
            double_moving_f = lyniar_double_moving_k(df, double_moving[0])
            print("The model is lyniar\n",double_moving_f[1],"\nThe forecast is made using ", mse[best_mse], "with k=", double_moving[0],
                  "\nThe forecast is:\n",
                  double_moving_f[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- regresia":
            e_f_reg = e_reg(df)
            print("The model is exponntial\n",e_f_reg[1],"\nThe forecast is made using ", mse[best_mse], "\nThe forecast is:\n",
                  e_f_reg[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- double moving":
            e_f_m = e_double_moving_f_ft_k_mse(df)
            print("The model is exponntial\n",e_f_m[1],"\nThe forecast is made using ", mse[best_mse],"with k=",e_f_m[2], "\nThe forecast is:\n", e_f_m[0],
                  "\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- model brown":
            e_f_b = e_brown_f_ft_alpha_mse(df)
            print("The model is exponntial\n",e_f_b[1],"\nThe forecast is made using ",mse[best_mse],"with alpha=",e_f_b[2],"\nThe forecast is:\n",e_f_b[0],
                  "\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- model holt":
            e_f_h = e_holt_f_ft_alpha_beta_mse(df)
            print("The model is exponntial\n",e_f_h[1],"\nThe forecast is made using ",mse[best_mse],"with alpha=",e_f_h[2],"and beta=",e_f_h[3],
                  "\nThe forecast is:\n",e_f_h[0],"\nThe mse is:", best_mse)
    else:
        c1 = c(df)
        season_avg = decomposition(c1)
        season_exponntial = season_exponntial_alpha_mse(season_avg)
        season_moving = season_moving_avg_k_mse(season_avg)
        season_fixed_avg = season_avg_fixed_mse(season_avg)
        season_double = season_lyniar_double_moving_k_mse(season_avg)
        season_holt = season_lyniar_holt_alpha_beta_mse(season_avg)
        season_brown = season_lyniar_brown_alpha_mse(season_avg)
        season_reg = season_lyniar_reg_mse(season_avg)
        season_square = season_square_mse(season_avg)
        k = best_k_mse(df)
        alpha = best_alpha_mse(df)
        regular_avg = fixed_avg_mse(df)
        reg_mse = lyniar_reg_mse(df)
        holt = best_alpha_beta_mse(df)
        print("A few more seconds and you have the forecast!!\n")
        brown = brown_best_alpha_mse(df)
        double_moving = double_best_k_mse(df)
        square = square_mse(df)
        reg_e = e_reg_mse(df)
        e_double = e_double_moving_f_ft_k_mse(df)
        e_brown = e_brown_f_ft_alpha_mse(df)
        e_holt = e_holt_f_ft_alpha_beta_mse(df)
        mse[season_square] = "season-squer-regresia"
        mse[season_double[1]] = "season- double moving avg"
        mse[season_holt[2]] = "season- model holt"
        mse[season_brown[1]] = "season- model brown"
        mse[season_reg] = "season-lyniar-regresia"
        mse[season_exponntial[1]] = "season- exponntial smooting"
        mse[season_moving[1]] = "season- moving average"
        mse[season_fixed_avg] = "season- regular average"
        mse[reg_e] = "exponntial- regresia"
        mse[e_double[3]] = "exponntial- double moving"
        mse[e_brown[3]] = "exponntial- model brown"
        mse[e_holt[4]] = "exponntial- model holt"
        mse[square] = "square regresia"
        mse[reg_mse] = "rgresia"
        mse[holt[2]] = "model holt"
        mse[brown[1]] = "model brown"
        mse[double_moving[1]] = "double moving avg"
        mse[regular_avg] = "regular avg"
        mse[k[1]] = "mooving avg"
        mse[alpha[1]] = "exponntial smooting"
        best_mse = min(mse)
        ab = reg1(season_avg)
        if (mse[best_mse] == "season- regular average" or mse[
            best_mse] == "season- moving average" or mse[
            best_mse] == "season- exponntial smooting" or float(ab[1]) == 0):
            print("The model is season- fixed\n",f(season_avg, c1)[1],"\nThe forecast is:\n", f(season_avg, c1)[0],"\nseasonality coefficients:\n")
            for i in range(1,13):
                print(i,":",c1[i-1],"\n")
            print("The mse is:", best_mse)
        elif mse[best_mse] == "season- double moving avg" or mse[
            best_mse] == "season- model holt" or mse[best_mse] == "season- model brown" or \
                mse[best_mse] == "season-lyniar-regresia":
            print("The model is season- lyniar\n",f1(season_avg, c1)[1],"\nThe forecast is:\n", f1(season_avg, c1)[0],"\nseasonality coefficients:\n")
            for i in range(1,13):
                print(i,":",c1[i-1],"\n")
            print("The mse is:", best_mse)
        elif (mse[best_mse] == "regular avg") or (mse[best_mse] == "mooving avg" and k[1] == regular_avg):
            regular_f = fixed_avg(df)
            print("The model is fixed\n",regular_f[1],"\nThe forecast is made using regular avg","\nThe forecast is:\n",regular_f[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "mooving avg":
            moving_f = fixed_moving_k(df, k[0])
            print("The model is fixed\n",moving_f[1],"\nThe forecast is made using ", mse[best_mse], "with k=", k[0],
                  "\nThe forecast is:\n", moving_f[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial smooting":
            smooting_f = fixed_exponntial_smooting(df, alpha[0])
            print("The model is fixed\n",smooting_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", alpha[0],
                  "\nThe forecast is:\n", smooting_f[0], "\nThe mse is:", best_mse)
        elif mse[best_mse] == "rgresia":
            reg1_f = lyniar_reg(df)
            print("The model is lyniar\n",reg1_f[1],"\nThe forecast is made using ",mse[best_mse],"\nThe forecast is:\n",reg1_f[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "model holt":
            holt_f = lyniar_exponntial_smooting(df, holt[0], holt[1])
            print("The model is lyniar\n",holt_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=", holt[0],
                  " and beta=", holt[1],"\nThe forecast is:\n",
                  holt_f[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "model brown":
            brown_f = lyniar_double_exponntial_smooting(df, brown[0])
            print("The model is lyniar\n", brown_f[1],"\nThe forecast is made using ", mse[best_mse], "with alpha=",brown[0],
                  "\nThe forecast is:\n",brown_f[0],"\n The mse is:", best_mse)
        elif mse[best_mse] == "double moving avg":
            double_moving_f = lyniar_double_moving_k(df, double_moving[0])
            print("The model is lyniar\n",double_moving_f[1],"\nThe forecast is made using ", mse[best_mse], "with k=",double_moving[0],
                  "\nThe forecast is:\n",double_moving_f[0],"\nThe mse is:",best_mse)
        elif mse[best_mse] == "square regresia":
            square_F = square_f(df)
            print("The model is square\n",square_F[1],"\nThe forecast is made using ",mse[best_mse],"\nThe forecast is:\n",square_F[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- regresia":
            e_f_reg = e_reg(df)
            print("The model is exponntial\n", e_f_reg[1],"\nThe forecast is made using ", mse[best_mse], "\nThe forecast is:\n",e_f_reg[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- double moving":
            e_f_m = e_double_moving_f_ft_k_mse(df)
            print("The model is exponntial\n",e_f_m[1],"\nThe forecast is made using ", mse[best_mse], "\nThe forecast is:\n",e_f_m[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- model brown":
            e_f_b = e_brown_f_ft_alpha_mse(df)
            print("The model is exponntial\n",e_f_b[1],"\nThe forecast is made using ",mse[best_mse],"\nThe forecast is:\n",e_f_b[0],"\nThe mse is:", best_mse)
        elif mse[best_mse] == "exponntial- model holt":
            e_f_h = e_holt_f_ft_alpha_beta_mse(df)
            print("The model is exponntial\n",e_f_h[1],"\nThe forecast is made using ",mse[best_mse],"\nThe forecast is:\n",e_f_h[0],"\nThe mse is:", best_mse)
        else:
            print("The model is season- squaer\n",f2(season_avg, c1)[1],"\nThe forecast is:\n",f2(season_avg, c1)[0],"\nseasonality coefficients:\n")
            for i in range(1,13):
                print(i,":",c1[i-1],"\n")
            print("The mse is:", best_mse)


o(df)



