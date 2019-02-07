using System;
using UnityEngine;



public class CFMtonemapperConfig: ScriptableObject
{
    #region runtime data collection
    [Serializable]
    public struct RuntimeData
    {
        /// curve parameters for curve function (ax^2+bx) / (cx^2+dx +1)
        /// Notice: original aces curve has an "e" in it, we eliminate it with division.
        public Vector3 aces_a;
        public Vector3 aces_b;
        public Vector3 aces_c;
        public Vector3 aces_d;

        // runtime white point this point
        public Vector3 whitePoint;

        // original ACES parameters:
        public const float ACES_A = 17.9285f;
        public const float ACES_B = 0.214285f;
        public const float ACES_C = 17.3571f;
        public const float ACES_D = 4.21428f;
        public const float ACES_WHITE_POINT = 50.0f;

        public RuntimeData(Vector3 _aces_a, Vector3 _aces_b, Vector3 _aces_c, Vector3 _aces_d, Vector3 _whitePoint)
        {
            aces_a = _aces_a;
            aces_b = _aces_b;
            aces_c = _aces_c;
            aces_d = _aces_d;
            whitePoint = _whitePoint;
        }

        /// <summary>
        /// Return the original ACES parameters:
        /// </summary>
        /// <returns></returns>
        public static RuntimeData Default()
        {
            return new RuntimeData( new Vector3(ACES_A, ACES_A, ACES_A),
                                    new Vector3(ACES_B, ACES_B, ACES_B),
                                    new Vector3(ACES_C, ACES_C, ACES_C),
                                    new Vector3(ACES_D, ACES_D, ACES_D),
                                    new Vector3(ACES_WHITE_POINT, ACES_WHITE_POINT, ACES_WHITE_POINT)
                );
        }
    }
    #endregion

    #region config data collection
    [Serializable]
    public class TonemappingSetting
    {
        public Vector2 toe = new Vector2(7.2f, 7.5f);
        public float shoulderAngle = 60.9f;
        public float shoulderShoot = 18.4f;
        public float toeStrength = 70.74f;
        public float toeLift = 28.8f;
        public float shoulderStrength = 97.3f;
        public float shoulderLift = 80.0f;
        public Vector2 whitePoint = new Vector2(6.0f, 1.0f);
        public float gamma = 1.0f;

        public void ConstrainParameter()
        {
            // clamp toe point:
            toe.x = Mathf.Max(toe.x, 0.0f);
            toe.y = Mathf.Clamp(toe.y, 0.0f, 100.0f);
            toeStrength = Mathf.Clamp(toeStrength, 0.1f, 99.9f);
            toeLift = Mathf.Clamp(toeLift, 0.1f, 99.9f);

            // clamp shoulder angle:
            shoulderAngle = Mathf.Clamp(shoulderAngle, 0.1f, 89.0f);
            float rad = Mathf.Deg2Rad * shoulderAngle;
            Vector2 shoulderDir = new Vector2(Mathf.Cos(rad), Mathf.Sin(rad));
            float slope = shoulderDir.y / shoulderDir.x;

            // constrain toe point:
            float toeSlope = toe.y / toe.x;
            if (slope < toeSlope)
            {
                // make the toe slope less equal to linear part slope, this is good for fitting result.
                toe.y = slope * toe.x;
            }
            // constrain shoulder shoot:
            shoulderShoot = Mathf.Max(shoulderShoot, 0.0f);
            shoulderStrength = Mathf.Clamp(shoulderStrength, 0.1f, 99.9f);
            shoulderLift = Mathf.Clamp(shoulderLift, 0.1f, 99.9f);
            whitePoint.x = Mathf.Clamp(whitePoint.x, 0.1f, 50.9f);
            whitePoint.y = Mathf.Clamp01(whitePoint.y);
            gamma = Mathf.Max(gamma, .01f);
        }
    }
    #endregion

    [SerializeField]
    public TonemappingSetting config;

    [SerializeField]
    public RuntimeData runtimeData = RuntimeData.Default();
    
    public float fitStandardDeviation { get { return m_fitStdDev; } }
    protected float m_fitStdDev = 0.0f;

    protected static TunableACEScurve m_tmpOriginalCurve = new TunableACEScurve();
    protected static ACEStonemappingCurve m_tmpAcesCurve = new ACEStonemappingCurve();
    protected static ACESfitting.Config fitConfig = ACESfitting.Config.Default();
}
