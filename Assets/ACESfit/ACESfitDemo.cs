using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class ACESfitDemo : MonoBehaviour {
    #region debug
    public bool drawEditorCurve = true;
    public bool drawRuntimeCurve = true;
    public bool drawImportantPoint = true;
    public bool drawFitSamples = true;
    public bool draw1OrderDrev = false;
    public bool draw2OrderDrev = false;
    public bool debugFit = false;
    public bool debugResetCurve = false;
    protected float pointGizmosSize = 0.02f;
    protected float sampleGizmosSize = 0.01f;
    protected Color sampleColor = new Color(1, 1, 1);
    protected Color _1orderDrevCol = new Color(1, .5f, 0);
    protected Color _2orderDrevCol = new Color(1, 1, 0);
    protected Color _editorCurveColor = new Color(1, 0, 0);
    protected Color _editorCurvePoint = new Color(1, 0, 0, 0.25f);

    protected Color handleColor = new Color(0, 0, 0);
    #endregion

    [SerializeField]
    public TunableACEScurve curve = new TunableACEScurve();
    public Vector2 toe = new Vector2(.05f, .05f);
    public float shootAngle = 45.0f;
    public float shootDistance = 10.0f;
    public Vector2 pWhite = new Vector2(5, 1);
    public Vector2 toeCtr = new Vector2(50.0f, 50.0f);
    public Vector2 shoulderCtr = new Vector2(50.0f, 50.0f);
    public float gamma = 1.0f;
    public float fitWhitePoint = 0.0f;

    [SerializeField]
    public ACESfitting.Config fitConfig = ACESfitting.Config.Default();

    public ACEStonemappingCurve acesCurve = new ACEStonemappingCurve();

    public MonotonicCubicBezier mCurve = new MonotonicCubicBezier();


    List<double> _tmpX;
    List<double> _tmpY;

    CurveDataSet data = null;

    const float A = 17.928571f;
    const float B = 0.214285f;
    const float C = 17.3571f;
    const float D = 4.21428f;
    const float E = 1.0f;
    const float whitePoint = 50.0f;

    void Start()
    {
        Shader.SetGlobalVector("param", new Vector4(A,B,C,D));
        Shader.SetGlobalFloat("whitePoint", whitePoint);
    }

    void Update () {
        curve.Toe = toe/100.0f;
        curve.ShootAngle = shootAngle;
        curve.Shoot = shootDistance/100.0f;
        curve.WhitePoint = pWhite;
        curve.ToeStrength = toeCtr.x / 100.0f;
        curve.ToeLift = toeCtr.y / 100.0f;
        curve.ShoulderStrength = shoulderCtr.x / 100.0f;
        curve.ShoulderLift = shoulderCtr.y / 100.0f;
        curve.Gamma = gamma;

        if (Input.GetKeyUp(KeyCode.Space) || debugFit)
        {
            debugFit = false;
            fitConfig.fitRange.Set(0, pWhite.x);
            ACESfitting.Fit(curve, ref acesCurve, fitConfig, out data, out _tmpX, out _tmpY);
            fitWhitePoint = acesCurve.CalculateWhitePoint();
            //Debug.Log(ACEStonemappingCurve.IsParamValid(acesCurve.a, acesCurve.b, acesCurve.c, acesCurve.d, fitConfig.fitRange.y, fitConfig.fitRange.x));
            Shader.SetGlobalVector("param", new Vector4(acesCurve.a, acesCurve.b, acesCurve.c, acesCurve.d));
            Shader.SetGlobalFloat("whitePoint", fitWhitePoint);
        }

        if (debugResetCurve)
        {
            debugResetCurve = false;
            acesCurve.Reset();
        }
    }


    private void OnDrawGizmos()
    {
        // fitting samples:
        if(_tmpX!=null&&_tmpY!=null && drawFitSamples)
        {
            Gizmos.color = sampleColor;

            for (int i = 0; i < _tmpX.Count; ++i)
            {
                Gizmos.DrawSphere(new Vector3((float)_tmpX[i], (float)_tmpY[i], 0.0f),sampleGizmosSize);
            }
        }

        // LDR curve
        Gizmos.color = Color.white;
        Gizmos.DrawLine(new Vector3(0, 0, 0), new Vector3(1, 1, 0));

        if (drawEditorCurve)
        {
            // draw important point
            if (curve != null && drawImportantPoint)
            {
                // toe:
                Gizmos.color = Color.red;
                Gizmos.DrawSphere(new Vector3(curve.Toe.x, curve.Toe.y, 0), pointGizmosSize);
                // shoulder:
                Gizmos.color = _editorCurvePoint;
                Gizmos.DrawSphere(new Vector3(curve.Shoulder.x, curve.Shoulder.y, 0), pointGizmosSize);

                // toe angle Point:
                Gizmos.color = _editorCurvePoint;
                Vector3 angleToe = new Vector3(curve.ToeAnglePoint.x, curve.ToeAnglePoint.y, 0);
                Gizmos.DrawSphere(angleToe, sampleGizmosSize);
                Gizmos.DrawLine(angleToe, curve.Toe);
                Gizmos.DrawLine(Vector3.zero, curve.Toe);
                //handle:
                Gizmos.color = handleColor;
                Vector3 toeStartHandle = Vector3.Lerp(Vector3.zero, angleToe, curve.ToeStrength);
                Gizmos.DrawLine(Vector3.zero, toeStartHandle);
                Gizmos.DrawSphere(toeStartHandle, sampleGizmosSize);
                Vector3 toeEndHandle = Vector3.Lerp(curve.Toe, angleToe, curve.ToeLift);
                Gizmos.DrawLine(curve.Toe, toeEndHandle);
                Gizmos.DrawSphere(toeEndHandle, sampleGizmosSize);

                Gizmos.color = _editorCurvePoint;
                // shoulder angle Point:
                Vector3 shoulderAngle = new Vector3(curve.ShoulderAnglePoint.x, curve.ShoulderAnglePoint.y, 0);
                Vector3 whitePoint = new Vector3(curve.WhitePoint.x, curve.WhitePoint.y, 0.0f);
                Gizmos.DrawSphere(shoulderAngle, sampleGizmosSize);
                Gizmos.DrawLine(shoulderAngle, whitePoint);
                Gizmos.DrawLine(shoulderAngle, curve.Shoulder);
                //handle:
                Gizmos.color = handleColor;
                Vector3 shoulderStartHandle = Vector3.Lerp(curve.Shoulder, shoulderAngle, curve.ShoulderStrength);
                Vector3 shoulder = new Vector3(curve.Shoulder.x, curve.Shoulder.y, 0);
                Gizmos.DrawLine(shoulder, shoulderStartHandle);
                Gizmos.DrawSphere(shoulderStartHandle, sampleGizmosSize);
                Vector3 shoulderEndHandle = Vector3.Lerp(whitePoint, shoulderAngle, curve.ShoulderLift);
                Gizmos.DrawLine(whitePoint, shoulderEndHandle);
                Gizmos.DrawSphere(shoulderEndHandle, sampleGizmosSize);

                // white point:
                Gizmos.color = Color.white;
                Gizmos.DrawSphere(new Vector3(curve.WhitePoint.x, curve.WhitePoint.y, 0), 0.025f);
            }

            // draw original data.
            Gizmos.color = Color.red;
            int drawCount = Mathf.RoundToInt( (fitConfig.fitRange.y - fitConfig.fitRange.x)* (float)100);
            float step = fitConfig.fitRange.y / (float)(drawCount);
            for (int i = 0; i < drawCount; ++i)
            {
                // original:
                float x = step * (float)i;
                float x_ = step * (float)(i+1);
                float y = curve.Eval(x);
                float y_ = curve.Eval(x_);
                Gizmos.DrawLine(new Vector3(x, y , 0), new Vector3(x_, y_ , 0));
            }

            if (data != null && draw1OrderDrev)
            {
                for (int i= 0; i < data.dataCount-1; ++i)
                {
                    Gizmos.color = _1orderDrevCol;
                    Gizmos.DrawLine(new Vector3((float)data.x[i], (float)data.dxdy[i], 0), new Vector3((float)data.x[i + 1], (float)data.dxdy[i + 1], 0));
                }
            }

            if(data != null && draw2OrderDrev)
            {
                for (int i = 0; i < data.dataCount - 1; ++i)
                {
                    Gizmos.color = _2orderDrevCol;
                    Gizmos.DrawLine(new Vector3((float)data.x[i], (float)data.dx2dy2[i], 0), new Vector3((float)data.x[i + 1], (float)data.dx2dy2[i + 1], 0));
                }
            }
        }
        
        if(drawRuntimeCurve)
        {
            int count = Mathf.RoundToInt((fitConfig.fitRange.y - fitConfig.fitRange.x) * (float)200);
            float step = (fitConfig.fitRange.y - fitConfig.fitRange.x )/ (float)count;

            for (int j = 0; j < count; ++j)
            {
                float _x = (float)j * step;
                Gizmos.color = Color.green;
                float _oy = (float)acesCurve.Eval(_x);
                float _oy_ = (float)acesCurve.Eval(_x+step);
                Gizmos.DrawLine(new Vector3(_x, Mathf.Pow(_oy, 1), 0),
                                 new Vector3(_x + step, Mathf.Pow(_oy_, 1), 0));
            }
        }
    }
    
}
